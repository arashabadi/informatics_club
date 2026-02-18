import React, { useEffect, useMemo, useRef, useState } from 'react';
import { useFrame } from '@react-three/fiber';
import { Line, Sphere, Html } from '@react-three/drei';
import * as THREE from 'three';
import { CellData, GameStage, KernelParams, KernelType } from '../types';
import { getTransitionProbs } from '../utils/simulation';

interface CellsProps {
  data: CellData[];
  activeCell: number | null;
  stage: GameStage;
  onCellClick: (id: number) => void;
  walkPath: number[];
  kernel: KernelType;
  kernelParams: KernelParams;
}

export const CellCloud: React.FC<CellsProps> = ({ data, activeCell, stage, onCellClick, walkPath, kernel, kernelParams }) => {
  const meshRef = useRef<THREE.InstancedMesh>(null);
  const tempObj = useMemo(() => new THREE.Object3D(), []);
  const [hoveredCell, setHoveredCell] = useState<number | null>(null);
  const walkPathSet = useMemo(() => new Set(walkPath), [walkPath]);
  const activeNeighborSet = useMemo(() => {
    if (activeCell === null) return null;
    const active = data[activeCell];
    if (!active) return null;
    return new Set(active.neighbors);
  }, [activeCell, data]);

  // Update instance positions and colors
  useFrame(() => {
    if (!meshRef.current) return;
    
    data.forEach((cell, i) => {
      tempObj.position.copy(cell.position);
      
      // Scale logic
      let scale = 0.15;
      if (activeCell === cell.id) scale = 0.4;
      if (walkPathSet.has(cell.id)) scale = 0.25;
      
      tempObj.scale.setScalar(scale);
      tempObj.updateMatrix();
      meshRef.current!.setMatrixAt(i, tempObj.matrix);

      // Color logic
      const color = new THREE.Color(cell.color);
      
      // Gray out others in later stages to focus
      const isNetworkStage = stage === GameStage.KNN_GRAPH || stage === GameStage.KERNEL_BIAS || stage === GameStage.TRANSITION_MATRIX;
      
      if (isNetworkStage && activeCell !== null && activeNeighborSet) {
         if (activeCell !== cell.id && !activeNeighborSet.has(cell.id)) {
            color.lerp(new THREE.Color('#1e293b'), 0.8); // dim
         }
      }
      
      if (walkPathSet.has(cell.id)) {
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
      const instanceId = e.instanceId;
      if (instanceId !== undefined) {
        setHoveredCell(instanceId);
      }
  };

  const handlePointerMove = (e: any) => {
      const instanceId = e.instanceId;
      if (instanceId !== undefined) {
        setHoveredCell(instanceId);
      }
  };
  
  const handlePointerOut = () => {
      document.body.style.cursor = 'default';
      setHoveredCell(null);
  };

  useEffect(() => {
    return () => {
      document.body.style.cursor = 'default';
    };
  }, []);

  const isNetworkStage = stage === GameStage.KNN_GRAPH || stage === GameStage.KERNEL_BIAS || stage === GameStage.TRANSITION_MATRIX;

  return (
    <group>
      <instancedMesh
        ref={meshRef}
        args={[undefined, undefined, data.length]}
        onClick={handleClick}
        onPointerOver={handlePointerOver}
        onPointerMove={handlePointerMove}
        onPointerOut={handlePointerOut}
      >
        <sphereGeometry args={[1, 16, 16]} />
        <meshStandardMaterial />
      </instancedMesh>

      {hoveredCell !== null && hoveredCell !== activeCell && data[hoveredCell] && (
        <Html
          position={[
            data[hoveredCell].position.x,
            data[hoveredCell].position.y + 0.45,
            data[hoveredCell].position.z,
          ]}
          center
          distanceFactor={12}
        >
          <div className="pointer-events-none rounded bg-slate-950/90 border border-slate-600 px-2 py-1 text-[10px] font-semibold text-slate-100 whitespace-nowrap">
            Cell #{hoveredCell}
          </div>
        </Html>
      )}

      {activeCell !== null && data[activeCell] && (
        <Html
          position={[
            data[activeCell].position.x,
            data[activeCell].position.y + 0.72,
            data[activeCell].position.z,
          ]}
          center
          distanceFactor={11}
        >
          <div className="pointer-events-none rounded bg-emerald-950/95 border border-emerald-400/70 px-2.5 py-1 text-[10px] font-bold text-emerald-200 whitespace-nowrap shadow-[0_0_10px_rgba(16,185,129,0.35)]">
            Selected Cell #{activeCell}
          </div>
        </Html>
      )}

      {/* Visualize Neighborhood Connections for Active Cell */}
      {isNetworkStage && activeCell !== null && (
        <Connections 
          data={data} 
          activeId={activeCell} 
          showDirection={stage === GameStage.KERNEL_BIAS || stage === GameStage.TRANSITION_MATRIX} 
          kernel={kernel}
          kernelParams={kernelParams}
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
const Connections: React.FC<{ data: CellData[], activeId: number, showDirection: boolean, kernel: KernelType, kernelParams: KernelParams }> = ({ data, activeId, showDirection, kernel, kernelParams }) => {
    const activeCell = data[activeId];
    const transitionProbs = useMemo(
      () => (activeCell ? getTransitionProbs(activeId, data, kernel, kernelParams) : []),
      [activeId, data, kernel, kernelParams, activeCell]
    );
    const maxProb = transitionProbs[0]?.prob ?? 1;
    const probByTarget = useMemo(() => {
      const probMap = new Map<number, number>();
      transitionProbs.forEach((entry) => probMap.set(entry.targetId, entry.prob));
      return probMap;
    }, [transitionProbs]);
    if (!activeCell) return null;

    return (
        <group>
            {activeCell.neighbors.map(nid => {
                const neighbor = data[nid];
                if (!neighbor) return null;
                const probability = probByTarget.get(nid) ?? 0;
                const normalized = maxProb > 0 ? probability / maxProb : 0;
                
                // Defaults
                let opacity = 0.12;
                let color = 'white';
                let thickness = 0.8;

                if (showDirection) {
                    if (kernel === 'Pseudotime') {
                        const dt = neighbor.pseudotime - activeCell.pseudotime;
                        if (dt > 0) {
                            opacity = 0.25 + normalized * 0.75;
                            color = '#4ade80';
                            thickness = 1 + normalized * 4;
                        } else {
                            opacity = 0.04;
                            color = '#64748b';
                            thickness = 0.4;
                        }
                    } 
                    else if (kernel === 'Velocity') {
                        const disp = new THREE.Vector3().subVectors(neighbor.position, activeCell.position).normalize();
                        const cos = activeCell.velocity.dot(disp);
                        if (cos > 0) {
                            opacity = Math.max(0.2, 0.3 + cos * normalized);
                            color = '#f472b6';
                            thickness = 1 + normalized * 4;
                        } else {
                            opacity = 0.04;
                            thickness = 0.5;
                        }
                    }
                    else if (kernel === 'CytoTRACE') {
                        const diff = activeCell.potency - neighbor.potency;
                        if (diff > 0) {
                             opacity = 0.2 + normalized * 0.8;
                             color = '#facc15';
                             thickness = 1 + normalized * 3.6;
                        } else {
                            opacity = 0.06;
                            color = '#713f12';
                            thickness = 0.5;
                        }
                    }
                    else if (kernel === 'RealTime') {
                         const dt = neighbor.pseudotime - activeCell.pseudotime;
                         const target = kernelParams.otTimeTarget;
                         const width = Math.max(0.03, target + 0.06);
                         const temporalMatch = Math.max(0, 1 - Math.abs(dt - target) / width);
                         if (temporalMatch > 0.15 && normalized > 0.25) {
                             opacity = 0.25 + normalized * temporalMatch * 0.75;
                             color = '#06b6d4';
                             thickness = 0.8 + normalized * temporalMatch * 5;
                         } else {
                             opacity = 0.02;
                             thickness = 0.3;
                         }
                    }
                    else if (kernel === 'Combined') {
                        const dt = neighbor.pseudotime - activeCell.pseudotime;
                        const disp = new THREE.Vector3().subVectors(neighbor.position, activeCell.position).normalize();
                        const cos = activeCell.velocity.dot(disp);
                        const directionBlend = ((cos + 1) / 2) * (dt > 0 ? 1 : 0.45);
                        opacity = 0.12 + normalized * directionBlend * 0.88;
                        color = '#818cf8';
                        thickness = 0.8 + normalized * 4;
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
    const points = useMemo(
      () => path.map(id => data[id]?.position).filter(Boolean) as THREE.Vector3[],
      [path, data]
    );
    if (points.length < 2) return null;
    return <Line points={points} color="white" lineWidth={3} transparent opacity={0.8} />;
};

export const TerminalHighlights: React.FC<{ data: CellData[], ids: number[] }> = ({ data, ids }) => {
    return (
        <group>
            {ids.map(id => {
                const cell = data[id];
                if (!cell) return null;
                return (
                <group key={id} position={cell.position}>
                     <Sphere args={[0.6, 16, 16]}>
                         <meshStandardMaterial color="#ef4444" emissive="#ef4444" emissiveIntensity={2} />
                     </Sphere>
                     <Html position={[0, 1, 0]}>
                         <div className="bg-red-600 text-white px-2 py-1 rounded text-xs font-bold whitespace-nowrap">
                             Terminal State
                         </div>
                     </Html>
                </group>
                );
            })}
        </group>
    )
}
