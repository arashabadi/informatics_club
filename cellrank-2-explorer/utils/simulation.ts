import * as THREE from 'three';
import { CellData, KernelParams, KernelType } from '../types';

// Constants
const NUM_CELLS = 1200;
const NEIGHBORS_K = 15; // Increased slightly for better connectivity in sparse kernels

export const DEFAULT_KERNEL_PARAMS: KernelParams = {
  pseudotimeBias: 10,
  velocitySigma: 0.35,
  cytoScale: 10,
  otEpsilon: 1,
  otTimeTarget: 0.05,
  alpha: 0.5,
};

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
  kernel: KernelType,
  params?: Partial<KernelParams>
): { targetId: number; prob: number }[] => {
  const resolved: KernelParams = {
    ...DEFAULT_KERNEL_PARAMS,
    ...params,
  };
  resolved.velocitySigma = Math.max(0.05, resolved.velocitySigma);
  resolved.otEpsilon = Math.max(0.05, resolved.otEpsilon);
  resolved.alpha = Math.min(1, Math.max(0, resolved.alpha));

  const current = cells[currentId];
  if (!current) return [];

  const candidates = current.neighbors.map(nid => cells[nid]);
  let weights: number[] = [];

  const computePseudotimeWeight = (c: CellData, n: CellData) => {
      const dt = n.pseudotime - c.pseudotime;
      return Math.max(1e-6, Math.exp(resolved.pseudotimeBias * dt));
  };

  const computeVelocityWeight = (c: CellData, n: CellData) => {
      const displacement = new THREE.Vector3().subVectors(n.position, c.position).normalize();
      const cosine = c.velocity.dot(displacement);
      return Math.max(1e-6, Math.exp(cosine / resolved.velocitySigma));
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
          return 1 / (1 + Math.exp(-resolved.cytoScale * diff));
      });
  }
  else if (kernel === 'RealTime') {
      // Optimal Transport approximation
      // Prefer neighbors in the "next" time bin (pseudotime + epsilon)
      // Penalize large distance
      weights = candidates.map(n => {
          const dt = n.pseudotime - current.pseudotime;
          // We want dt to be around otTimeTarget
          const timeError = Math.abs(dt - resolved.otTimeTarget);
          const dist = current.position.distanceTo(n.position);

          const cost = ((dist * dist) + (timeError * 12)) / resolved.otEpsilon;
          return Math.max(1e-6, Math.exp(-cost));
      });
  }
  else if (kernel === 'Combined') {
      // 50% Velocity + 50% Pseudotime
      weights = candidates.map(n => {
          const wVel = computeVelocityWeight(current, n);
          const wPseudo = computePseudotimeWeight(current, n);
          return resolved.alpha * wVel + (1 - resolved.alpha) * wPseudo;
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
  kernel: KernelType,
  params?: Partial<KernelParams>
): number => {
  const probs = getTransitionProbs(currentId, cells, kernel, params);
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

export const computeKernelCompareScores = (
  currentId: number,
  cells: CellData[],
  params?: Partial<KernelParams>
): { kernel: KernelType; score: number }[] => {
  const kernels: KernelType[] = ['Pseudotime', 'Velocity', 'CytoTRACE', 'RealTime', 'Combined'];
  return kernels.map((kernel) => {
    const probs = getTransitionProbs(currentId, cells, kernel, params);
    return {
      kernel,
      score: probs[0]?.prob ?? 0,
    };
  });
};

export const getInitialStates = (cells: CellData[]): number[] => {
    if (cells.length === 0) return [];
    const stemLike = cells
      .filter(c => c.type === 'Stem' || c.type === 'Progenitor')
      .sort((a, b) => a.pseudotime - b.pseudotime);

    const chosen = (stemLike.length > 0 ? stemLike : [...cells].sort((a, b) => a.pseudotime - b.pseudotime))
      .slice(0, 2)
      .map(c => c.id);

    return Array.from(new Set(chosen));
};

// Find macrostates (fake implementation for visual)
export const getMacrostates = (cells: CellData[]): number[] => {
    if (cells.length === 0) return [];
    const typeA = cells.filter(c => c.type === 'TypeA').sort((a,b) => b.pseudotime - a.pseudotime)[0] ?? cells[0];
    const typeB = cells.filter(c => c.type === 'TypeB').sort((a,b) => b.pseudotime - a.pseudotime)[0] ?? cells[cells.length - 1];
    if (typeA.id === typeB.id && cells.length > 1) {
      return [typeA.id, cells[cells.length - 1].id];
    }
    return [typeA.id, typeB.id];
}
