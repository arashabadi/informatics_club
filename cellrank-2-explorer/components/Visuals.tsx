import React, { useMemo, useRef, useState } from 'react';
import { useFrame } from '@react-three/fiber';
import { Line, Sphere, Html } from '@react-three/drei';
import * as THREE from 'three';
import { CellData, GameStage, KernelType } from '../types';

interface CellsProps {
  data: CellData[];
  activeCell: number | null;
  stage: GameStage;
  onCellClick: (id: number) => void;
  walkPath: number[];
  kernel: KernelType;
}

export const CellCloud: React.FC<CellsProps> = ({ data, activeCell, stage, onCellClick, walkPath, kernel }) => {
  const meshRef = useRef<THREE.InstancedMesh>(null);
  const tempObj = new THREE.Object3D();
  const hoverRef = useRef<number | null>(null);

  // Update instance positions and colors
  useFrame(() => {
    if (!meshRef.current) return;
    
    data.forEach((cell, i) => {
      tempObj.position.copy(cell.position);
      
      // Scale logic
      let scale = 0.15;
      if (activeCell === cell.id) scale = 0.4;
      if (walkPath.includes(cell.id)) scale = 0.25;
      
      tempObj.scale.setScalar(scale);
      tempObj.updateMatrix();
      meshRef.current!.setMatrixAt(i, tempObj.matrix);

      // Color logic
      const color = new THREE.Color(cell.color);
      
      // Gray out others in later stages to focus
      const isNetworkStage = stage === GameStage.KNN_GRAPH || stage === GameStage.KERNEL_BIAS || stage === GameStage.TRANSITION_MATRIX;
      
      if (isNetworkStage && activeCell !== null) {
         if (activeCell !== cell.id && !data[activeCell].neighbors.includes(cell.id)) {
            color.lerp(new THREE.Color('#1e293b'), 0.8); // dim
         }
      }
      
      if (walkPath.includes(cell.id)) {
          color.set('#ffffff'); // Path is white
          if (cell.id === walkPath[walkPath.length-1]) color.set('#10b981'); // Head is green
      }

      meshRef.current!.setColorAt(i, color);
    });
    meshRef.current.instanceMatrix.needsUpdate = true;
    if (meshRef.current.instanceColor) meshRef.current.instanceColor.needsUpdate = true;
  });

  // Handle raycasting for clicks
  const handleClick = (e: any) => {
    e.stopPropagation();
    const instanceId = e.instanceId;
    if (instanceId !== undefined) {
        onCellClick(instanceId);
    }
  };

  const handlePointerOver = (e: any) => {
      e.stopPropagation();
      document.body.style.cursor = 'pointer';
      hoverRef.current = e.instanceId;
  };
  
  const handlePointerOut = () => {
      document.body.style.cursor = 'default';
      hoverRef.current = null;
  };

  const isNetworkStage = stage === GameStage.KNN_GRAPH || stage === GameStage.KERNEL_BIAS || stage === GameStage.TRANSITION_MATRIX;

  return (
    <group>
      <instancedMesh
        ref={meshRef}
        args={[undefined, undefined, data.length]}
        onClick={handleClick}
        onPointerOver={handlePointerOver}
        onPointerOut={handlePointerOut}
      >
        <sphereGeometry args={[1, 16, 16]} />
        <meshStandardMaterial />
      </instancedMesh>

      {/* Visualize Neighborhood Connections for Active Cell */}
      {isNetworkStage && activeCell !== null && (
        <Connections 
          data={data} 
          activeId={activeCell} 
          showDirection={stage === GameStage.KERNEL_BIAS || stage === GameStage.TRANSITION_MATRIX} 
          kernel={kernel}
        />
      )}

      {/* Walk Path Line */}
      {walkPath.length > 1 && (
         <WalkLine data={data} path={walkPath} />
      )}
    </group>
  );
};

// Draws lines/arrows between active cell and neighbors
const Connections: React.FC<{ data: CellData[], activeId: number, showDirection: boolean, kernel: KernelType }> = ({ data, activeId, showDirection, kernel }) => {
    const activeCell = data[activeId];
    if (!activeCell) return null;

    return (
        <group>
            {activeCell.neighbors.map(nid => {
                const neighbor = data[nid];
                
                // Defaults
                let opacity = 0.3;
                let color = 'white';
                let thickness = 1;

                if (showDirection) {
                    if (kernel === 'Pseudotime') {
                        if (neighbor.pseudotime > activeCell.pseudotime) {
                            opacity = 1; color = '#4ade80'; thickness = 3; // Green forward
                        } else {
                            opacity = 0.1; color = '#94a3b8';
                        }
                    } 
                    else if (kernel === 'Velocity') {
                        const disp = new THREE.Vector3().subVectors(neighbor.position, activeCell.position).normalize();
                        const cos = activeCell.velocity.dot(disp);
                        if (cos > 0) {
                            opacity = 0.5 + cos * 0.5; color = '#f472b6'; thickness = 2 + cos * 3; // Pink velocity
                        } else {
                            opacity = 0.1;
                        }
                    }
                    else if (kernel === 'CytoTRACE') {
                        // High potency -> Low potency
                        const diff = activeCell.potency - neighbor.potency;
                        if (diff > 0) {
                             opacity = 1; color = '#facc15'; thickness = 3; // Yellow flow
                        } else {
                            opacity = 0.1; color = '#94a3b8';
                        }
                    }
                    else if (kernel === 'RealTime') {
                         // Similar to Pseudotime but tighter
                         const dt = neighbor.pseudotime - activeCell.pseudotime;
                         if (dt > 0 && dt < 0.1) {
                             opacity = 1; color = '#06b6d4'; thickness = 3; // Cyan tight temporal coupling
                         } else {
                             opacity = 0.1;
                         }
                    }
                    else if (kernel === 'Combined') {
                        // Mix of Velocity (Pink) and Pseudotime (Green) -> Purple/White?
                        // Let's use a unique Indigo
                         if (neighbor.pseudotime > activeCell.pseudotime) {
                             opacity = 1; color = '#818cf8'; thickness = 3;
                         } else {
                             opacity = 0.1;
                         }
                    }
                }

                return (
                    <Line
                        key={nid}
                        points={[activeCell.position, neighbor.position]}
                        color={color}
                        transparent
                        opacity={opacity}
                        lineWidth={thickness}
                    />
                )
            })}
            
            {/* Show velocity vector for Velocity Kernel */}
            {showDirection && (kernel === 'Velocity' || kernel === 'Combined') && (
                <arrowHelper 
                    args={[
                        activeCell.velocity, 
                        activeCell.position, 
                        2.5, 
                        0xff00ff
                    ]} 
                />
            )}
        </group>
    );
};

// Draws the random walk path
const WalkLine: React.FC<{ data: CellData[], path: number[] }> = ({ data, path }) => {
    const points = useMemo(() => path.map(id => data[id].position), [path, data]);
    return <Line points={points} color="white" lineWidth={3} transparent opacity={0.8} />;
};

export const TerminalHighlights: React.FC<{ data: CellData[], ids: number[] }> = ({ data, ids }) => {
    return (
        <group>
            {ids.map(id => (
                <group key={id} position={data[id].position}>
                     <Sphere args={[0.6, 16, 16]}>
                         <meshStandardMaterial color="#ef4444" emissive="#ef4444" emissiveIntensity={2} />
                     </Sphere>
                     <Html position={[0, 1, 0]}>
                         <div className="bg-red-600 text-white px-2 py-1 rounded text-xs font-bold whitespace-nowrap">
                             Terminal State
                         </div>
                     </Html>
                </group>
            ))}
        </group>
    )
}