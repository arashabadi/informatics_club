import React, { useState, useEffect } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Stars, Environment } from '@react-three/drei';
import { CellData, GameStage, KernelType } from '../types';
import { generateManifold, nextStep, getMacrostates } from '../utils/simulation';
import { CellCloud, TerminalHighlights } from './Visuals';
import { Interface } from './Interface';

export const Explorer3D: React.FC<{ onBack: () => void }> = ({ onBack }) => {
  const [stage, setStage] = useState<GameStage>(GameStage.INTRO);
  const [cells, setCells] = useState<CellData[]>([]);
  const [activeCell, setActiveCell] = useState<number | null>(null);
  const [walkPath, setWalkPath] = useState<number[]>([]);
  const [kernel, setKernel] = useState<KernelType>('Pseudotime');
  const [isWalking, setIsWalking] = useState(false);
  const [terminalStates, setTerminalStates] = useState<number[]>([]);

  useEffect(() => {
    const data = generateManifold();
    setCells(data);
    setTerminalStates(getMacrostates(data));
    setActiveCell(10); 
  }, []);

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'ArrowRight') setStage(s => Math.min(s + 1, GameStage.MACROSTATES));
      else if (e.key === 'ArrowLeft') setStage(s => Math.max(s - 1, GameStage.INTRO));
      else if (e.key === ' ') if (stage === GameStage.RANDOM_WALK) handleWalkStep();
    };
    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [stage, activeCell, cells, kernel]);

  useEffect(() => {
    setWalkPath([]);
    setIsWalking(false);
    if (stage === GameStage.RANDOM_WALK && activeCell === null) setActiveCell(10); 
  }, [stage, kernel]);

  const handleCellClick = (id: number) => {
    if (stage === GameStage.MACROSTATES) return;
    setActiveCell(id);
    setWalkPath([id]);
  };

  const handleWalkStep = () => {
    if (activeCell === null) return;
    let currentPath = walkPath.length > 0 ? [...walkPath] : [activeCell];
    const currentId = currentPath[currentPath.length - 1];
    const nextId = nextStep(currentId, cells, kernel);
    currentPath.push(nextId);
    setWalkPath(currentPath);
  };

  const handleResetWalk = () => {
    if (activeCell !== null) setWalkPath([activeCell]);
    else setWalkPath([]);
  };

  return (
    <div className="w-full h-screen bg-slate-900 relative">
      <button 
        onClick={onBack}
        className="absolute top-4 left-4 z-50 bg-slate-800 text-white px-3 py-1 rounded border border-slate-600 hover:bg-slate-700"
      >
        Exit to Menu
      </button>

      <Canvas gl={{ preserveDrawingBuffer: true }} camera={{ position: [0, 0, 18], fov: 45 }}>
        <color attach="background" args={['#0f172a']} />
        <ambientLight intensity={0.5} />
        <pointLight position={[10, 10, 10]} intensity={1} />
        <Stars radius={100} depth={50} count={5000} factor={4} saturation={0} fade speed={1} />
        <Environment preset="city" />
        <group rotation={[0, -Math.PI / 6, 0]}>
            <CellCloud 
                data={cells} 
                activeCell={activeCell} 
                stage={stage} 
                onCellClick={handleCellClick}
                walkPath={walkPath}
                kernel={kernel}
            />
            {stage === GameStage.MACROSTATES && <TerminalHighlights data={cells} ids={terminalStates} />}
        </group>
        <OrbitControls enablePan={true} enableZoom={true} enableRotate={true} autoRotate={stage === GameStage.INTRO} autoRotateSpeed={0.5} />
      </Canvas>

      <Interface 
        stage={stage}
        setStage={setStage}
        kernel={kernel}
        setKernel={setKernel}
        onWalk={handleWalkStep}
        onResetWalk={handleResetWalk}
        isWalking={isWalking}
        cells={cells}
        activeCell={activeCell}
      />
    </div>
  );
};
