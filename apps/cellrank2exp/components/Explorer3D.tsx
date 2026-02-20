import React, { useState, useEffect, useCallback } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Stars, AdaptiveDpr } from '@react-three/drei';
import * as THREE from 'three';
import { CellData, GameStage, KernelParams, KernelType } from '../types';
import { DEFAULT_KERNEL_PARAMS, generateManifold, nextStep, getInitialStates, getMacrostates } from '../utils/simulation';
import { CellCloud, InitialHighlights, TerminalHighlights } from './Visuals';
import { Interface } from './Interface';

const MAX_WALK_LENGTH = 800;

export const Explorer3D: React.FC<{ onBack: () => void }> = ({ onBack }) => {
  const [stage, setStage] = useState<GameStage>(GameStage.INTRO);
  const [cells, setCells] = useState<CellData[]>([]);
  const [activeCell, setActiveCell] = useState<number | null>(null);
  const [walkPath, setWalkPath] = useState<number[]>([]);
  const [kernel, setKernel] = useState<KernelType>('Pseudotime');
  const [kernelParams, setKernelParams] = useState<KernelParams>({ ...DEFAULT_KERNEL_PARAMS });
  const [isWalking, setIsWalking] = useState(false);
  const [initialStates, setInitialStates] = useState<number[]>([]);
  const [terminalStates, setTerminalStates] = useState<number[]>([]);
  const [matrixFocus, setMatrixFocus] = useState<{ sourceId: number; targetId: number; value: number } | null>(null);
  const [navigationMode, setNavigationMode] = useState<'ROTATE' | 'PAN'>('ROTATE');
  const [isInitializing, setIsInitializing] = useState(true);
  const [initError, setInitError] = useState<string | null>(null);
  const [initVersion, setInitVersion] = useState(0);

  useEffect(() => {
    let isCancelled = false;
    let frameId = 0;
    setIsInitializing(true);
    setInitError(null);

    frameId = window.requestAnimationFrame(() => {
      try {
        const data = generateManifold();
        const initialIds = getInitialStates(data);
        const terminalIds = getMacrostates(data);
        if (isCancelled) return;
        setCells(data);
        setInitialStates(initialIds);
        setTerminalStates(terminalIds);
        setActiveCell(initialIds[0] ?? (data.length > 10 ? 10 : data[0]?.id ?? null));
      } catch (error) {
        if (isCancelled) return;
        const message = error instanceof Error ? error.message : 'Unknown explorer initialization error';
        setInitError(message);
        setCells([]);
        setInitialStates([]);
        setTerminalStates([]);
        setActiveCell(null);
      } finally {
        if (!isCancelled) setIsInitializing(false);
      }
    });

    return () => {
      isCancelled = true;
      window.cancelAnimationFrame(frameId);
    };
  }, [initVersion]);

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'ArrowRight') setStage(s => Math.min(s + 1, GameStage.MACROSTATES));
      else if (e.key === 'ArrowLeft') setStage(s => Math.max(s - 1, GameStage.INTRO));
      else if (e.key === ' ' && stage === GameStage.RANDOM_WALK) {
        if (activeCell === null || cells.length === 0) return;
        setWalkPath((prevPath) => {
          const currentPath = prevPath.length > 0 ? [...prevPath] : [activeCell];
          const currentId = currentPath[currentPath.length - 1];
          const nextId = nextStep(currentId, cells, kernel, kernelParams);
          currentPath.push(nextId);
          if (currentPath.length > MAX_WALK_LENGTH) {
            currentPath.splice(0, currentPath.length - MAX_WALK_LENGTH);
          }
          return currentPath;
        });
      }
    };
    window.addEventListener('keydown', handleKeyDown);
    return () => window.removeEventListener('keydown', handleKeyDown);
  }, [stage, activeCell, cells, kernel, kernelParams]);

  useEffect(() => {
    setWalkPath([]);
    setIsWalking(false);
    if (stage === GameStage.RANDOM_WALK && activeCell === null && cells.length > 0) {
      setActiveCell(Math.min(10, cells.length - 1));
    }
    if (stage !== GameStage.TRANSITION_MATRIX) {
      setMatrixFocus(null);
    }
  }, [stage, kernel, activeCell, cells.length]);

  const handleCellClick = (id: number) => {
    if (id < 0 || id >= cells.length) return;
    if (stage === GameStage.MACROSTATES) return;
    setMatrixFocus(null);
    setActiveCell(id);
    setWalkPath([id]);
  };

  const handleMatrixCellSelect = useCallback((sourceId: number, targetId: number, value: number) => {
    if (sourceId < 0 || sourceId >= cells.length) return;
    if (targetId < 0 || targetId >= cells.length) return;
    setMatrixFocus({ sourceId, targetId, value });
    setActiveCell(sourceId);
    setWalkPath([]);
  }, [cells.length]);

  const handleWalkStep = useCallback(() => {
    if (activeCell === null || cells.length === 0) return;
    setWalkPath((prevPath) => {
      const currentPath = prevPath.length > 0 ? [...prevPath] : [activeCell];
      const currentId = currentPath[currentPath.length - 1];
      const nextId = nextStep(currentId, cells, kernel, kernelParams);
      currentPath.push(nextId);
      if (currentPath.length > MAX_WALK_LENGTH) {
        currentPath.splice(0, currentPath.length - MAX_WALK_LENGTH);
      }
      return currentPath;
    });
  }, [activeCell, cells, kernel, kernelParams]);

  const handleToggleWalk = () => {
    if (activeCell === null) return;
    setIsWalking((prev) => !prev);
  };

  const handleWalkBurst = (steps: number) => {
    if (steps <= 0) return;
    for (let i = 0; i < steps; i += 1) {
      handleWalkStep();
    }
  };

  const handleWalkBack = () => {
    setWalkPath((prevPath) => {
      if (prevPath.length <= 1) return prevPath;
      return prevPath.slice(0, -1);
    });
  };

  const handleKernelParamChange = (key: keyof KernelParams, value: number) => {
    setKernelParams((prev) => ({
      ...prev,
      [key]: value,
    }));
  };

  useEffect(() => {
    if (!isWalking || stage !== GameStage.RANDOM_WALK) return;
    const timer = window.setInterval(() => {
      handleWalkStep();
    }, 350);
    return () => window.clearInterval(timer);
  }, [isWalking, stage, handleWalkStep]);

  const handleResetWalk = () => {
    if (activeCell !== null) setWalkPath([activeCell]);
    else setWalkPath([]);
  };

  useEffect(() => {
    if (activeCell === null) return;
    if (stage !== GameStage.RANDOM_WALK) return;
    setWalkPath([activeCell]);
  }, [kernelParams, activeCell, stage]);

  return (
    <div className="w-full h-screen bg-slate-900 relative">
      {isInitializing && (
        <div className="absolute inset-0 z-40 flex items-center justify-center bg-slate-950/75 text-slate-200 text-sm">
          Initializing Explorer 3D…
        </div>
      )}

      {initError && (
        <div className="absolute top-16 left-4 right-4 z-50 rounded border border-red-500/40 bg-red-950/85 px-3 py-2 text-sm text-red-200">
          <div>Explorer failed to initialize: {initError}</div>
          <button
            onClick={() => setInitVersion((prev) => prev + 1)}
            className="mt-2 rounded border border-red-300/50 bg-red-900/60 px-2 py-1 text-xs text-red-100 hover:bg-red-800/70"
          >
            Retry Explorer Init
          </button>
        </div>
      )}

      <Canvas
        dpr={[1, 1.5]}
        gl={{ preserveDrawingBuffer: true, antialias: true, alpha: false }}
        performance={{ min: 0.6 }}
        camera={{ position: [0, 0, 18], fov: 45 }}
        onCreated={({ gl }) => {
          gl.setClearColor('#020617');
        }}
      >
        <color attach="background" args={['#0f172a']} />
        <AdaptiveDpr pixelated />
        <ambientLight intensity={0.75} />
        <hemisphereLight args={['#93c5fd', '#0f172a', 0.45]} />
        <pointLight position={[10, 10, 10]} intensity={1.2} />
        <Stars radius={90} depth={40} count={1800} factor={3} saturation={0} fade speed={0.7} />
        <gridHelper args={[28, 28, '#334155', '#1e293b']} position={[0, -4.5, 0]} />
        <group rotation={[0, -Math.PI / 6, 0]}>
            <CellCloud 
                data={cells} 
                activeCell={activeCell} 
                matrixFocus={matrixFocus}
                stage={stage} 
                onCellClick={handleCellClick}
                walkPath={walkPath}
                kernel={kernel}
                kernelParams={kernelParams}
            />
            {(stage === GameStage.INTRO || stage === GameStage.MACROSTATES) && <InitialHighlights data={cells} ids={initialStates} />}
            {stage === GameStage.MACROSTATES && <TerminalHighlights data={cells} ids={terminalStates} />}
        </group>
        <OrbitControls
          enablePan={true}
          enableZoom={true}
          enableRotate={true}
          screenSpacePanning={true}
          panSpeed={1.8}
          zoomSpeed={1.2}
          rotateSpeed={0.7}
          enableDamping={true}
          dampingFactor={0.08}
          mouseButtons={{
            LEFT: navigationMode === 'PAN' ? THREE.MOUSE.PAN : THREE.MOUSE.ROTATE,
            MIDDLE: THREE.MOUSE.DOLLY,
            RIGHT: navigationMode === 'PAN' ? THREE.MOUSE.ROTATE : THREE.MOUSE.PAN,
          }}
          autoRotate={stage === GameStage.INTRO}
          autoRotateSpeed={0.5}
        />
      </Canvas>

      <Interface 
        onExit={onBack}
        initialStates={initialStates}
        terminalStates={terminalStates}
        matrixFocus={matrixFocus}
        navigationMode={navigationMode}
        setNavigationMode={setNavigationMode}
        stage={stage}
        setStage={setStage}
        kernel={kernel}
        setKernel={setKernel}
        onWalk={handleWalkStep}
        onWalkBack={handleWalkBack}
        onResetWalk={handleResetWalk}
        onMatrixCellSelect={handleMatrixCellSelect}
        onToggleWalk={handleToggleWalk}
        onWalkBurst={handleWalkBurst}
        isWalking={isWalking}
        cells={cells}
        activeCell={activeCell}
        walkPath={walkPath}
        kernelParams={kernelParams}
        onKernelParamChange={handleKernelParamChange}
      />
    </div>
  );
};
