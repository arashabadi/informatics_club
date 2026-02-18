import * as THREE from 'three';
import { CellData, KernelType } from '../types';

// Constants
const NUM_CELLS = 1200;
const NEIGHBORS_K = 15; // Increased slightly for better connectivity in sparse kernels

// Helper to generate a Y-shaped bifurcation manifold
export const generateManifold = (): CellData[] => {
  const cells: CellData[] = [];
  
  for (let i = 0; i < NUM_CELLS; i++) {
    // 0 to 1 progress
    const t = Math.random();
    
    // Y-branching logic
    let x, y, z;
    let type: CellData['type'] = 'Stem';
    let color = '#e2e8f0'; // slate-200
    
    // Add some noise
    const noise = () => (Math.random() - 0.5) * 1.5;

    if (t < 0.4) {
      // Stem/Progenitor trunk
      x = t * 10;
      y = noise();
      z = noise();
      type = t < 0.1 ? 'Stem' : 'Progenitor';
      color = t < 0.1 ? '#fde047' : '#cbd5e1'; // Yellow start
    } else {
      // Branching
      const branch = Math.random() > 0.5 ? 1 : -1;
      const branchProgress = (t - 0.4) / 0.6; // 0 to 1 along branch
      
      x = 4 + branchProgress * 10;
      y = branch * (branchProgress * 8) + noise();
      z = noise();
      
      if (branch === 1) {
        type = 'TypeA';
        // Gradient to Red
        color = `rgb(${100 + branchProgress * 155}, 100, 100)`; 
      } else {
        type = 'TypeB';
        // Gradient to Blue
        color = `rgb(100, 100, ${100 + branchProgress * 155})`;
      }
    }

    // Velocity vector (approximate tangent)
    const velocity = new THREE.Vector3(1, 0, 0); 
    if (type === 'TypeA') velocity.set(1, 0.8, 0).normalize();
    if (type === 'TypeB') velocity.set(1, -0.8, 0).normalize();
    
    // Add some stochasticity to velocity
    velocity.x += (Math.random() - 0.5) * 0.2;
    velocity.y += (Math.random() - 0.5) * 0.2;
    velocity.normalize();

    // Potency: Inverse of pseudotime (Stem has high potency) + noise
    // In real CytoTRACE, this is gene count.
    const potency = (1 - t) + (Math.random() * 0.1); 

    cells.push({
      id: i,
      position: new THREE.Vector3(x - 5, y, z), // Center scene roughly
      color,
      pseudotime: t,
      potency,
      velocity,
      neighbors: [],
      type
    });
  }

  // Compute KNN (Naive O(N^2) for simplicity in browser, manageable for 1200 cells)
  cells.forEach(cell => {
    const distances = cells.map(other => ({
      id: other.id,
      dist: cell.position.distanceTo(other.position)
    }));
    
    // Sort by distance and pick top K (exclude self)
    distances.sort((a, b) => a.dist - b.dist);
    cell.neighbors = distances.slice(1, NEIGHBORS_K + 1).map(d => d.id);
  });

  return cells;
};

// Calculate transition probabilities for a specific cell
export const getTransitionProbs = (
  currentId: number,
  cells: CellData[],
  kernel: KernelType
): { targetId: number; prob: number }[] => {
  const current = cells[currentId];
  if (!current) return [];

  const candidates = current.neighbors.map(nid => cells[nid]);
  let weights: number[] = [];

  const computePseudotimeWeight = (c: CellData, n: CellData) => {
      const dt = n.pseudotime - c.pseudotime;
      if (dt > 0) return 1.0 + dt * 10; 
      return 0.05;
  };

  const computeVelocityWeight = (c: CellData, n: CellData) => {
      const displacement = new THREE.Vector3().subVectors(n.position, c.position).normalize();
      const cosine = c.velocity.dot(displacement);
      return Math.max(0.01, (cosine + 1) ** 2);
  };

  if (kernel === 'Pseudotime') {
     weights = candidates.map(n => computePseudotimeWeight(current, n));
  } 
  else if (kernel === 'Velocity') {
    weights = candidates.map(n => computeVelocityWeight(current, n));
  }
  else if (kernel === 'CytoTRACE') {
      // Flow from High Potency -> Low Potency
      // w_ij ~ logistic(potency_i - potency_j)
      weights = candidates.map(n => {
          const diff = current.potency - n.potency; // Positive if current is stem-like and neighbor is diff-like
          // Sigmoid-like scaling
          return 1 / (1 + Math.exp(-10 * diff)); 
      });
  }
  else if (kernel === 'RealTime') {
      // Optimal Transport approximation
      // Prefer neighbors in the "next" time bin (pseudotime + epsilon)
      // Penalize large distance
      const epsilon = 0.05; // Time step size
      weights = candidates.map(n => {
          const dt = n.pseudotime - current.pseudotime;
          // We want dt to be around epsilon
          const timeError = Math.abs(dt - epsilon);
          const dist = current.position.distanceTo(n.position);
          
          // Cost = distance^2 + lambda * time_deviation
          // Weight ~ exp(-Cost)
          const cost = (dist * 0.5) + (timeError * 20);
          return Math.exp(-cost);
      });
  }
  else if (kernel === 'Combined') {
      // 50% Velocity + 50% Pseudotime
      weights = candidates.map(n => {
          const wVel = computeVelocityWeight(current, n);
          const wPseudo = computePseudotimeWeight(current, n);
          // Simple additive combination (in reality usually linear comb of transition matrices)
          return 0.5 * wVel + 0.5 * wPseudo;
      });
  }

  // Normalize weights
  const totalWeight = weights.reduce((a, b) => a + b, 0);
  
  if (totalWeight === 0) return candidates.map(c => ({ targetId: c.id, prob: 1/candidates.length }));

  // Return objects
  return candidates.map((c, i) => ({
      targetId: c.id,
      prob: weights[i] / totalWeight
  })).sort((a, b) => b.prob - a.prob);
};

// Simulate a random walk step based on a kernel
export const nextStep = (
  currentId: number, 
  cells: CellData[], 
  kernel: KernelType
): number => {
  const probs = getTransitionProbs(currentId, cells, kernel);
  if (probs.length === 0) return currentId;

  // Sample
  const r = Math.random();
  let cumSum = 0;
  for (const p of probs) {
    cumSum += p.prob;
    if (r <= cumSum) return p.targetId;
  }
  
  return probs[probs.length - 1].targetId;
};

// Find macrostates (fake implementation for visual)
export const getMacrostates = (cells: CellData[]): number[] => {
    const typeA = cells.filter(c => c.type === 'TypeA').sort((a,b) => b.pseudotime - a.pseudotime)[0];
    const typeB = cells.filter(c => c.type === 'TypeB').sort((a,b) => b.pseudotime - a.pseudotime)[0];
    return [typeA.id, typeB.id];
}