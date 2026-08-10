import React, { useEffect, useMemo, useState } from 'react';
import { CellData, GameStage, KernelParams, KernelType } from '../types';
import { ArrowRight, Play, Pause, FastForward, RefreshCcw, Info, GitGraph, Activity, Zap, Grid, Camera, Loader2, ChevronDown, LogOut, Undo2 } from 'lucide-react';
import { computeKernelCompareScores, getTransitionProbs } from '../utils/simulation';
// @ts-ignore
import html2canvas from 'html2canvas';

interface InterfaceProps {
  onExit: () => void;
  initialStates: number[];
  terminalStates: number[];
  matrixFocus: { sourceId: number; targetId: number; value: number } | null;
  navigationMode: 'ROTATE' | 'PAN';
  setNavigationMode: (mode: 'ROTATE' | 'PAN') => void;
  stage: GameStage;
  setStage: (s: GameStage) => void;
  kernel: KernelType;
  setKernel: (k: KernelType) => void;
  onWalk: () => void;
  onWalkBack: () => void;
  onResetWalk: () => void;
  onMatrixCellSelect: (sourceId: number, targetId: number, value: number) => void;
  onToggleWalk: () => void;
  onWalkBurst: (steps: number) => void;
  isWalking: boolean;
  cells: CellData[];
  activeCell: number | null;
  walkPath: number[];
  kernelParams: KernelParams;
  onKernelParamChange: (key: keyof KernelParams, value: number) => void;
}

const Formula = ({ children }: { children?: React.ReactNode }) => (
    <div className="font-mono bg-slate-800 p-2 rounded text-emerald-400 text-sm my-2 border border-slate-700 overflow-x-auto">
        {children}
    </div>
);

// Helper for math text
const M = ({ children }: { children?: React.ReactNode }) => (
    <span className="font-serif italic">{children}</span>
);

const ParamSlider = ({
  label,
  value,
  min,
  max,
  step,
  onChange,
}: {
  label: string;
  value: number;
  min: number;
  max: number;
  step: number;
  onChange: (value: number) => void;
}) => (
  <label className="block">
    <div className="flex justify-between items-center text-[11px] text-slate-400 mb-1">
      <span>{label}</span>
      <span className="text-slate-200 font-mono">{value.toFixed(step < 1 ? 2 : 1)}</span>
    </div>
    <input
      className="w-full accent-emerald-500"
      type="range"
      min={min}
      max={max}
      step={step}
      value={value}
      onChange={(event) => onChange(parseFloat(event.target.value))}
    />
  </label>
);

export const Interface: React.FC<InterfaceProps> = ({ 
  onExit,
  initialStates,
  terminalStates,
  matrixFocus,
  navigationMode,
  setNavigationMode,
  stage,
  setStage,
  kernel,
  setKernel,
  onWalk,
  onWalkBack,
  onResetWalk,
  onMatrixCellSelect,
  onToggleWalk,
  onWalkBurst,
  isWalking,
  cells,
  activeCell,
  walkPath,
  kernelParams,
  onKernelParamChange,
}) => {
  const [isCapturing, setIsCapturing] = useState(false);
  const [showKernelMenu, setShowKernelMenu] = useState(false);
  const [hoveredMatrixEntry, setHoveredMatrixEntry] = useState<{ source: number; target: number; value: number } | null>(null);
  const stageLabels = ['Manifold', 'KNN Graph', 'Kernel', 'Matrix', 'Walk', 'Macrostates'];
  const selectedCell = activeCell !== null ? cells[activeCell] : null;
  const selectedCellRole = selectedCell
    ? `${initialStates.includes(selectedCell.id) ? 'Initial' : ''}${initialStates.includes(selectedCell.id) && terminalStates.includes(selectedCell.id) ? ' + ' : ''}${terminalStates.includes(selectedCell.id) ? 'Terminal' : ''}` || 'Intermediate'
    : null;
  const uniqueVisited = walkPath.length > 0 ? new Set(walkPath).size : 0;
  const activeProbs = useMemo(
    () => (activeCell !== null ? getTransitionProbs(activeCell, cells, kernel, kernelParams) : []),
    [activeCell, cells, kernel, kernelParams]
  );
  const kernelScores = activeCell !== null ? computeKernelCompareScores(activeCell, cells, kernelParams) : [];
  const totalScore = kernelScores.reduce((acc, score) => acc + score.score, 0);
  const normalizedScores = kernelScores.map((entry) => ({
    ...entry,
    normalized: totalScore > 0 ? entry.score / totalScore : 0,
  }));
  const matrixPreview = useMemo(() => {
    if (activeCell === null || activeProbs.length === 0) return null;
    const indices = Array.from(new Set([activeCell, ...activeProbs.slice(0, 11).map((entry) => entry.targetId)])).slice(0, 12);
    if (indices.length < 2) return null;

    const rows = indices.map((sourceId) => {
      const rowProbs = getTransitionProbs(sourceId, cells, kernel, kernelParams);
      const probMap = new Map<number, number>(rowProbs.map((entry) => [entry.targetId, entry.prob]));
      return indices.map((targetId) => probMap.get(targetId) ?? 0);
    });

    const maxValue = rows.flat().reduce((max, value) => Math.max(max, value), 0);
    return { indices, rows, maxValue: maxValue > 0 ? maxValue : 1 };
  }, [activeCell, activeProbs, cells, kernel, kernelParams]);

  const matrixPythonSnippet = useMemo(() => {
    const kernelInit: Record<KernelType, string> = {
      Pseudotime: 'kernel = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")',
      Velocity: 'kernel = cr.kernels.VelocityKernel(adata)',
      CytoTRACE: 'kernel = cr.kernels.CytoTRACEKernel(adata, potency_key="ct_potency")',
      RealTime: 'kernel = cr.kernels.RealTimeKernel(adata, time_key="time")',
      Combined: 'kernel = 0.6 * vk + 0.4 * pk',
    };

    return `import cellrank as cr
import matplotlib.pyplot as plt
import numpy as np

${kernelInit[kernel]}
kernel.compute_transition_matrix()

T = kernel.transition_matrix  # scipy.sparse.csr_matrix
print("shape:", T.shape, "nnz:", T.nnz)

idx = np.random.choice(T.shape[0], 40, replace=False)
block = T[idx][:, idx].A

plt.figure(figsize=(4, 4))
plt.imshow(block, cmap="viridis")
plt.colorbar(label="P(i→j)")
plt.title("Transition matrix block")`;
  }, [kernel]);

  useEffect(() => {
    setHoveredMatrixEntry(null);
  }, [stage, activeCell, kernel]);

  const renderKernelControls = () => {
    if (stage !== GameStage.KERNEL_BIAS) return null;

    return (
      <div className="space-y-2 mt-3 border border-slate-700 bg-slate-800/70 rounded-lg p-3">
        <div className="text-[11px] uppercase tracking-wider text-slate-400">Live Formula Controls</div>
        {(kernel === 'Pseudotime' || kernel === 'Combined') && (
          <ParamSlider
            label="Pseudotime Bias b"
            value={kernelParams.pseudotimeBias}
            min={1}
            max={20}
            step={0.5}
            onChange={(value) => onKernelParamChange('pseudotimeBias', value)}
          />
        )}
        {(kernel === 'Velocity' || kernel === 'Combined') && (
          <ParamSlider
            label="Velocity Sigma σ"
            value={kernelParams.velocitySigma}
            min={0.1}
            max={2}
            step={0.05}
            onChange={(value) => onKernelParamChange('velocitySigma', value)}
          />
        )}
        {kernel === 'CytoTRACE' && (
          <ParamSlider
            label="CytoTRACE Scale"
            value={kernelParams.cytoScale}
            min={1}
            max={20}
            step={0.5}
            onChange={(value) => onKernelParamChange('cytoScale', value)}
          />
        )}
        {kernel === 'RealTime' && (
          <>
            <ParamSlider
              label="OT Epsilon ε"
              value={kernelParams.otEpsilon}
              min={0.1}
              max={3}
              step={0.1}
              onChange={(value) => onKernelParamChange('otEpsilon', value)}
            />
            <ParamSlider
              label="Target Δt"
              value={kernelParams.otTimeTarget}
              min={0}
              max={0.25}
              step={0.01}
              onChange={(value) => onKernelParamChange('otTimeTarget', value)}
            />
          </>
        )}
        {kernel === 'Combined' && (
          <ParamSlider
            label="Combination α (Velocity)"
            value={kernelParams.alpha}
            min={0}
            max={1}
            step={0.05}
            onChange={(value) => onKernelParamChange('alpha', value)}
          />
        )}
      </div>
    );
  };

  const renderKernelCompare = () => {
    if (activeCell === null || (stage !== GameStage.KERNEL_BIAS && stage !== GameStage.TRANSITION_MATRIX)) {
      return null;
    }

    return (
      <div className="space-y-2 mt-3 border border-slate-700 bg-slate-800/70 rounded-lg p-3">
        <div className="text-[11px] uppercase tracking-wider text-slate-400">Kernel Quick Compare (Cell #{activeCell})</div>
        {normalizedScores.map((entry) => (
          <div key={entry.kernel}>
            <div className="flex justify-between text-[11px] mb-1">
              <span className="text-slate-300">{entry.kernel}</span>
              <span className="text-slate-200 font-mono">{(entry.normalized * 100).toFixed(1)}%</span>
            </div>
            <div className="h-2 rounded bg-slate-700 overflow-hidden">
              <div className="h-full bg-gradient-to-r from-cyan-500 to-blue-500" style={{ width: `${entry.normalized * 100}%` }} />
            </div>
          </div>
        ))}
      </div>
    );
  };

  const handleSnapshot = async () => {
      setIsCapturing(true);
      try {
          const canvas = await html2canvas(document.body, {
              backgroundColor: '#0f172a',
              ignoreElements: (element: any) => element.id === 'ui-controls'
          });
          
          const link = document.createElement('a');
          link.download = `cellrank-slide-${stage}-${Date.now()}.png`;
          link.href = canvas.toDataURL('image/png');
          link.click();
      } catch (e) {
          console.error("Snapshot failed", e);
      }
      setIsCapturing(false);
  };
  
  const kernels: { id: KernelType; label: string; color: string }[] = [
      { id: 'Pseudotime', label: 'Pseudotime', color: 'blue' },
      { id: 'Velocity', label: 'Velocity', color: 'purple' },
      { id: 'CytoTRACE', label: 'CytoTRACE', color: 'yellow' },
      { id: 'RealTime', label: 'RealTime (OT)', color: 'cyan' },
      { id: 'Combined', label: 'Combined (Vel+Ps)', color: 'indigo' },
  ];

  const renderContent = () => {
    switch (stage) {
      case GameStage.INTRO:
        return (
          <div className="space-y-4">
            <h2 className="text-2xl font-bold text-blue-400 flex items-center gap-2"><GitGraph /> The Manifold</h2>
            <p className="text-slate-300">
              In single-cell biology, we often view cells as points in a high-dimensional space. 
              Here, we see a low-dimensional projection (like UMAP) of a differentiating lineage.
            </p>
            <p className="text-slate-300">
              The <strong>Yellow</strong> cells are Stem cells. They differentiate into two distinct lineages: 
              <strong> Red (Type A)</strong> and <strong>Blue (Type B)</strong>.
            </p>
            <div className="rounded border border-emerald-500/40 bg-emerald-950/30 p-2 text-xs text-emerald-200">
              Initial states are shown as green markers (source-like, low pseudotime). IDs: {initialStates.length > 0 ? initialStates.map((id) => `#${id}`).join(', ') : 'N/A'}.
            </div>
          </div>
        );
      case GameStage.KNN_GRAPH:
        return (
          <div className="space-y-4">
            <h2 className="text-2xl font-bold text-blue-400 flex items-center gap-2"><Activity /> The KNN Graph</h2>
            <p className="text-slate-300">
              CellRank models cell dynamics as a graph. For every cell <M>i</M>, we find its <M>k</M>-nearest neighbors based on transcriptomic similarity.
            </p>
            <Formula>
              <M>N(i)</M> = &#123; <M>j</M> | <M>j</M> is among <M>k</M>-nearest neighbors of <M>i</M> &#125;
            </Formula>
            <p className="text-sm text-yellow-200 bg-yellow-900/30 p-2 rounded">
              <Info className="inline w-4 h-4 mr-1"/> Interaction: Click any cell to see its connections. Notice connections are symmetric and undirected currently.
            </p>
          </div>
        );
      case GameStage.KERNEL_BIAS:
        return (
          <div className="space-y-4">
            <h2 className="text-2xl font-bold text-purple-400 flex items-center gap-2"><Zap /> The Kernel</h2>
            <p className="text-slate-300">
              Select a Kernel to see how it biases the graph edges.
            </p>
            
            <div className="relative">
                <button 
                    onClick={() => setShowKernelMenu(!showKernelMenu)}
                    className="w-full bg-slate-800 border border-slate-700 rounded-lg p-3 flex justify-between items-center text-sm font-bold hover:bg-slate-700 transition"
                >
                    <span>{kernels.find(k => k.id === kernel)?.label}</span>
                    <ChevronDown size={16} />
                </button>
                
                {showKernelMenu && (
                    <div className="absolute top-full left-0 w-full mt-2 bg-slate-900 border border-slate-700 rounded-lg shadow-xl z-50 overflow-hidden">
                        {kernels.map(k => (
                            <button
                                key={k.id}
                                onClick={() => { setKernel(k.id); setShowKernelMenu(false); }}
                                className={`w-full text-left px-4 py-2 text-sm hover:bg-slate-800 flex items-center gap-2 ${kernel === k.id ? 'text-' + k.color + '-400 font-bold' : 'text-slate-400'}`}
                            >
                                <span className={`w-2 h-2 rounded-full bg-${k.color}-500`}></span>
                                {k.label}
                            </button>
                        ))}
                    </div>
                )}
            </div>

            <div className="min-h-[100px]">
                {kernel === 'Pseudotime' && (
                    <>
                        <Formula><M>P(j|i)</M> ∝ exp(<M>b</M> · (<M>t<sub>j</sub></M> - <M>t<sub>i</sub></M>))</Formula>
                        <p className="text-xs text-slate-400">Green Edge: <M>j</M> is older (forward).</p>
                    </>
                )}
                {kernel === 'Velocity' && (
                    <>
                        <Formula><M>P(j|i)</M> ∝ CosSim(<M>v<sub>i</sub></M>, <M>x<sub>j</sub></M> - <M>x<sub>i</sub></M>)</Formula>
                         <p className="text-xs text-slate-400">Pink Edge: Aligns with velocity.</p>
                    </>
                )}
                {kernel === 'CytoTRACE' && (
                    <>
                        <Formula><M>P(j|i)</M> ∝ 1 / (1 + exp(-<M>potency<sub>i</sub></M> + <M>potency<sub>j</sub></M>))</Formula>
                        <p className="text-xs text-slate-400">Yellow Edge: Flow from High to Low Potency.</p>
                    </>
                )}
                {kernel === 'RealTime' && (
                     <>
                        <Formula><M>P(j|i)</M> ∝ exp(-Cost(<M>t</M>, <M>dist</M>))</Formula>
                        <p className="text-xs text-slate-400">Cyan Edge: Optimal Transport coupling.</p>
                    </>
                )}
                 {kernel === 'Combined' && (
                     <>
                        <Formula><M>P</M> = 0.5·<M>P<sub>vel</sub></M> + 0.5·<M>P<sub>pseudo</sub></M></Formula>
                        <p className="text-xs text-slate-400">Indigo Edge: Consensus of multiple kernels.</p>
                    </>
                )}
            </div>

             <p className="text-sm text-yellow-200 bg-yellow-900/30 p-2 rounded">
              <Info className="inline w-4 h-4 mr-1"/> Click a cell to see biased edges.
            </p>
          </div>
        );
      case GameStage.TRANSITION_MATRIX:
        const probs = activeProbs;
        
        return (
          <div className="space-y-4">
            <h2 className="text-2xl font-bold text-green-400 flex items-center gap-2"><Grid /> Transition Matrix</h2>
            <p className="text-slate-300">
              Row-stochastic <strong>Transition Matrix <M>T</M></strong> based on {kernel} kernel.
            </p>
            <Formula>
                <M>T<sub>ij</sub></M> = <M>P</M>(cell<M><sub>i</sub></M> &#8594; cell<M><sub>j</sub></M>)
            </Formula>
            
            {activeCell !== null && (
                <div className="bg-slate-800 p-3 rounded-lg border border-slate-600">
                    <h3 className="text-xs font-bold text-slate-400 uppercase mb-2">
                        Local Transition Row <M>T<sub>{activeCell}, :</sub></M>
                    </h3>
                    <div className="space-y-1 max-h-40 overflow-y-auto">
                        {probs.map((p, i) => (
                            <div key={i} className="flex items-center gap-2 text-xs">
                                <span className="w-8 text-right text-slate-500">#{p.targetId}</span>
                                <div className="flex-1 h-3 bg-slate-700 rounded-full overflow-hidden">
                                    <div 
                                        className="h-full bg-green-500" 
                                        style={{ width: `${p.prob * 100}%` }} 
                                    />
                                </div>
                                <span className="w-10 text-right text-green-300">{(p.prob * 100).toFixed(1)}%</span>
                            </div>
                        ))}
                    </div>
                </div>
            )}

            {matrixPreview && (
                <div className="bg-slate-800 p-3 rounded-lg border border-slate-600">
                    <h3 className="text-xs font-bold text-slate-400 uppercase mb-2">
                        Matrix Block View (Python-like)
                    </h3>
                    <p className="text-[11px] text-slate-300 mb-2">
                      Row index <span className="text-cyan-300 font-mono">i</span> = source cell where a walk starts.
                      Column index <span className="text-amber-300 font-mono">j</span> = target cell after one step.
                      Click any box to highlight both cells in 3D.
                    </p>
                     <div className="overflow-x-auto">
                        <div
                            className="inline-grid gap-[2px] rounded-md bg-slate-900/70 p-2 border border-slate-700"
                            style={{ gridTemplateColumns: `repeat(${matrixPreview.indices.length}, minmax(10px, 1fr))` }}
                        >
                            {matrixPreview.rows.map((row, rowIdx) =>
                                row.map((value, colIdx) => {
                                  const sourceId = matrixPreview.indices[rowIdx];
                                  const targetId = matrixPreview.indices[colIdx];
                                  const isFocused = matrixFocus?.sourceId === sourceId && matrixFocus?.targetId === targetId;
                                  return (
                                    <div
                                        key={`${rowIdx}-${colIdx}`}
                                        className={`w-3.5 h-3.5 rounded-[2px] cursor-crosshair ${isFocused ? 'ring-2 ring-cyan-300/90' : ''}`}
                                        style={{ backgroundColor: `rgba(34, 197, 94, ${Math.max(0.07, value / matrixPreview.maxValue)})` }}
                                        title={`T[${matrixPreview.indices[rowIdx]}, ${matrixPreview.indices[colIdx]}] = ${value.toFixed(4)}`}
                                        onMouseEnter={() => setHoveredMatrixEntry({
                                            source: sourceId,
                                            target: targetId,
                                            value,
                                        })}
                                        onMouseLeave={() => setHoveredMatrixEntry((prev) => (prev && prev.source === matrixPreview.indices[rowIdx] && prev.target === matrixPreview.indices[colIdx] ? null : prev))}
                                        onClick={() => onMatrixCellSelect(sourceId, targetId, value)}
                                        onKeyDown={(event) => {
                                          if (event.key === 'Enter' || event.key === ' ') {
                                            event.preventDefault();
                                            onMatrixCellSelect(sourceId, targetId, value);
                                          }
                                        }}
                                        role="button"
                                        tabIndex={0}
                                        aria-label={`Select matrix entry T[${sourceId}, ${targetId}]`}
                                    />
                                  );
                                })
                            )}
                        </div>
                    </div>
                    <div className="text-[11px] text-emerald-300 font-mono mt-2">
                        {hoveredMatrixEntry
                          ? `T[${hoveredMatrixEntry.source}, ${hoveredMatrixEntry.target}] = ${hoveredMatrixEntry.value.toFixed(4)}`
                          : 'Hover a box to inspect probability values'}
                    </div>
                    {matrixFocus && (
                      <div className="text-[11px] text-cyan-200 font-mono mt-1">
                        Focused pair: row/source #{matrixFocus.sourceId} → column/target #{matrixFocus.targetId} | T = {matrixFocus.value.toFixed(4)}
                      </div>
                    )}
                    <p className="text-[11px] text-slate-500 mt-2 font-mono break-all">
                        Cell order: [{matrixPreview.indices.join(', ')}]
                    </p>
                </div>
            )}

            <div className="bg-slate-950/80 p-3 rounded-lg border border-slate-700">
                <h3 className="text-xs font-bold text-slate-400 uppercase mb-2">
                    How this looks in Python
                </h3>
                <pre className="text-[11px] leading-relaxed text-slate-200 overflow-x-auto font-mono">
{matrixPythonSnippet}
                </pre>
            </div>
            
            <p className="text-sm text-yellow-200 bg-yellow-900/30 p-2 rounded mt-2">
                <Info className="inline w-4 h-4 mr-1"/> Click cells to inspect probabilities.
            </p>
          </div>
        );
      case GameStage.RANDOM_WALK:
        return (
           <div className="space-y-4">
            <h2 className="text-2xl font-bold text-orange-400 flex items-center gap-2"><Play /> Random Walks</h2>
            <p className="text-slate-300">
              Simulate a Random Walk using the <strong>{kernel} Kernel</strong>.
            </p>
            <div className="flex flex-wrap gap-2">
                <button 
                    onClick={onWalk}
                    disabled={isWalking}
                    className="bg-orange-600 hover:bg-orange-500 disabled:opacity-50 text-white px-4 py-2 rounded font-bold flex items-center gap-2"
                >
                    <Play size={16} /> Simulate Step
                </button>
                <button
                    onClick={onWalkBack}
                    disabled={isWalking || walkPath.length <= 1}
                    className="bg-amber-700 hover:bg-amber-600 disabled:opacity-50 text-white px-4 py-2 rounded font-bold flex items-center gap-2"
                >
                    <Undo2 size={16} /> Step Back
                </button>
                <button 
                    onClick={onToggleWalk}
                    className="bg-purple-600 hover:bg-purple-500 text-white px-4 py-2 rounded font-bold flex items-center gap-2"
                >
                    {isWalking ? <Pause size={16} /> : <Play size={16} />}
                    {isWalking ? 'Pause Auto' : 'Auto Walk'}
                </button>
                <button 
                    onClick={() => onWalkBurst(10)}
                    disabled={isWalking}
                    className="bg-sky-700 hover:bg-sky-600 disabled:opacity-50 text-white px-4 py-2 rounded font-bold flex items-center gap-2"
                >
                    <FastForward size={16} /> +10 Steps
                </button>
                 <button 
                    onClick={onResetWalk}
                    className="bg-slate-700 hover:bg-slate-600 text-white px-4 py-2 rounded font-bold flex items-center gap-2"
                >
                    <RefreshCcw size={16} /> Reset
                </button>
            </div>
            <div className="grid grid-cols-2 gap-2">
                <div className="bg-slate-800/80 border border-slate-700 rounded px-3 py-2 text-xs">
                    <div className="text-slate-400">Walk length</div>
                    <div className="text-orange-300 font-mono">{walkPath.length}</div>
                </div>
                <div className="bg-slate-800/80 border border-slate-700 rounded px-3 py-2 text-xs">
                    <div className="text-slate-400">Unique cells visited</div>
                    <div className="text-orange-300 font-mono">{uniqueVisited}</div>
                </div>
            </div>
             <p className="text-sm text-yellow-200 bg-yellow-900/30 p-2 rounded">
              <Info className="inline w-4 h-4 mr-1"/> Watch the particle follow the kernel's bias.
            </p>
           </div>
        );
       case GameStage.MACROSTATES:
        return (
            <div className="space-y-4">
            <h2 className="text-2xl font-bold text-red-500">Macrostates, Initial vs Terminal</h2>
            <p className="text-slate-300">
              CellRank groups cells into coarse dynamical states (macrostates). A practical view is:
              <strong> initial states</strong> act as sources, <strong>terminal states</strong> act as sinks.
            </p>
            <p className="text-slate-300">
              Terminal states are regions where random walks accumulate. Initial states are early cells that mostly flow outward toward multiple fates.
            </p>
             <div className="grid grid-cols-1 gap-2">
                <div className="p-3 border border-emerald-500/30 bg-emerald-900/20 rounded">
                    <h3 className="font-bold text-emerald-300 mb-1 text-sm">Initial States (Sources)</h3>
                    <p className="text-xs text-emerald-100">IDs: {initialStates.length > 0 ? initialStates.map((id) => `#${id}`).join(', ') : 'N/A'}</p>
                </div>
                <div className="p-3 border border-red-500/30 bg-red-900/20 rounded">
                    <h3 className="font-bold text-red-300 mb-1 text-sm">Terminal States (Sinks)</h3>
                    <p className="text-xs text-red-100">IDs: {terminalStates.length > 0 ? terminalStates.map((id) => `#${id}`).join(', ') : 'N/A'}</p>
                </div>
             </div>
             <div className="mt-1 p-4 border border-red-500/30 bg-red-900/20 rounded">
                <h3 className="font-bold text-red-400 mb-2">Fate Probability Concept</h3>
                <p className="text-sm">
                    <M>P</M>(fate<M><sub>A</sub></M> | cell<M><sub>i</sub></M>) is the probability that a random walk starting at <M>i</M> ends in Terminal State A.
                </p>
                <p className="text-xs text-slate-300 mt-2">
                    In practice, high terminal-state probabilities indicate commitment; mixed probabilities indicate transitional/plastic cells.
                </p>
            </div>
           </div>
        )
      default:
        return null;
    }
  };

  return (
    <div className="absolute inset-0 pointer-events-none flex flex-col">
      {/* Header */}
      <div className="bg-slate-900/90 backdrop-blur border-b border-slate-700 px-4 py-3 pointer-events-auto flex justify-between items-center">
        <div>
            <h1 className="text-xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-blue-400 to-emerald-400">
                CellRank 2 Explorer
            </h1>
            <p className="text-xs text-slate-400">Interactive Tutorial for Biologists</p>
        </div>
        <div className="flex items-center gap-3">
             <button
                onClick={onExit}
                className="inline-flex items-center gap-1.5 rounded-md border border-slate-600 bg-slate-800 px-2.5 py-1.5 text-xs text-slate-100 hover:bg-slate-700"
                title="Exit to Menu"
             >
                <LogOut className="w-3.5 h-3.5" />
                Exit
             </button>
             <div className="hidden md:flex items-center gap-1">
                <button
                  onClick={() => setNavigationMode('ROTATE')}
                  className={`rounded border px-2 py-1 text-[11px] ${navigationMode === 'ROTATE' ? 'border-blue-400/70 bg-blue-500/20 text-blue-100' : 'border-slate-700 bg-slate-800 text-slate-300'}`}
                  title="Left-drag rotates, right-drag pans"
                >
                  Rotate Mode
                </button>
                <button
                  onClick={() => setNavigationMode('PAN')}
                  className={`rounded border px-2 py-1 text-[11px] ${navigationMode === 'PAN' ? 'border-blue-400/70 bg-blue-500/20 text-blue-100' : 'border-slate-700 bg-slate-800 text-slate-300'}`}
                  title="Left-drag pans, right-drag rotates"
                >
                  Pan Mode
                </button>
             </div>
             <div className="hidden xl:flex text-xs text-slate-500 gap-2 items-center">
                <span className="border border-slate-700 px-1 rounded bg-slate-800">wheel/pinch</span>
                <span>zoom</span>
                <span className="border border-slate-700 px-1 rounded bg-slate-800">mode buttons</span>
                <span>choose drag behavior</span>
             </div>
             
             <button 
                onClick={handleSnapshot}
                disabled={isCapturing}
                className="text-slate-400 hover:text-white transition-colors"
                title="Save Snapshot for Slides"
             >
                {isCapturing ? <Loader2 className="animate-spin w-5 h-5"/> : <Camera className="w-5 h-5" />}
             </button>
             
            <div className="text-xs text-slate-500 font-mono">
                {stage + 1} / 6
            </div>
        </div>
      </div>

      <div className="px-4 md:px-6 pt-2 pb-1 pointer-events-auto">
        <div className="max-w-3xl flex flex-nowrap gap-2 overflow-x-auto pb-1">
            {stageLabels.map((label, idx) => (
                <button
                    key={label}
                    onClick={() => setStage(idx as GameStage)}
                    className={`shrink-0 text-[11px] px-2.5 py-1 rounded-full border transition ${
                        stage === idx
                            ? 'bg-blue-500/20 border-blue-400/60 text-blue-200'
                            : 'bg-slate-800/70 border-slate-700 text-slate-400 hover:text-slate-200'
                    }`}
                >
                    {idx + 1}. {label}
                </button>
            ))}
        </div>
      </div>

      {/* Main Content Card */}
      <div className="flex-1 min-h-0 px-4 md:px-6 pb-2 md:pb-4 pointer-events-none">
        <div
          className="pointer-events-auto max-w-md h-full max-h-full overflow-y-auto overscroll-y-contain pr-1"
          onWheel={(event) => event.stopPropagation()}
          onTouchMove={(event) => event.stopPropagation()}
        >
          <div className="bg-slate-900/95 backdrop-blur border border-slate-700 p-4 md:p-6 rounded-xl shadow-2xl">
            {renderContent()}
            {renderKernelControls()}
            {renderKernelCompare()}

            {selectedCell && (
                <div className="mt-4 border border-slate-700 bg-slate-800/70 rounded-lg p-3 text-xs grid grid-cols-2 gap-2">
                    <div>
                        <div className="text-slate-500">Selected Cell</div>
                        <div className="text-slate-100 font-mono">#{selectedCell.id}</div>
                    </div>
                    <div>
                        <div className="text-slate-500">Type</div>
                        <div className="text-slate-100">{selectedCell.type}</div>
                    </div>
                    <div>
                        <div className="text-slate-500">State Role</div>
                        <div className="text-cyan-300">{selectedCellRole}</div>
                    </div>
                    <div>
                        <div className="text-slate-500">Pseudotime</div>
                        <div className="text-emerald-300 font-mono">{selectedCell.pseudotime.toFixed(3)}</div>
                    </div>
                    <div>
                        <div className="text-slate-500">Potency</div>
                        <div className="text-yellow-300 font-mono">{selectedCell.potency.toFixed(3)}</div>
                    </div>
                </div>
            )}
            
            <div className="mt-6 flex justify-between border-t border-slate-800 pt-4" id="ui-controls">
                <button 
                    disabled={stage === 0}
                    onClick={() => setStage(stage - 1)}
                    className="text-slate-400 hover:text-white text-sm font-medium disabled:opacity-30"
                >
                    &larr; Back
                </button>
                <button 
                    disabled={stage === 5}
                    onClick={() => setStage(stage + 1)}
                    className="bg-blue-600 hover:bg-blue-500 text-white px-4 py-2 rounded-lg text-sm font-bold flex items-center gap-2 disabled:opacity-30 disabled:cursor-not-allowed transition-colors"
                >
                    Next <ArrowRight size={16} />
                </button>
            </div>
          </div>
        </div>
      </div>
      
      {/* Legend Footer */}
      <div className="hidden md:flex p-4 bg-slate-900/80 backdrop-blur text-xs gap-6 justify-center text-slate-400 border-t border-slate-800 pointer-events-auto">
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-yellow-400"></span> Stem</div>
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-slate-300"></span> Progenitor</div>
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-red-500"></span> Type A</div>
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-blue-500"></span> Type B</div>
      </div>
    </div>
  );
};
