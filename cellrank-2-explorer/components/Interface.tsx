import React, { useState } from 'react';
import { CellData, GameStage, KernelType } from '../types';
import { ArrowRight, Play, Pause, FastForward, RefreshCcw, Info, GitGraph, Activity, Zap, Grid, Camera, Loader2, ChevronDown } from 'lucide-react';
import { getTransitionProbs } from '../utils/simulation';
// @ts-ignore
import html2canvas from 'html2canvas';

interface InterfaceProps {
  stage: GameStage;
  setStage: (s: GameStage) => void;
  kernel: KernelType;
  setKernel: (k: KernelType) => void;
  onWalk: () => void;
  onResetWalk: () => void;
  onToggleWalk: () => void;
  onWalkBurst: (steps: number) => void;
  isWalking: boolean;
  cells: CellData[];
  activeCell: number | null;
  walkPath: number[];
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

export const Interface: React.FC<InterfaceProps> = ({ 
  stage, setStage, kernel, setKernel, onWalk, onResetWalk, onToggleWalk, onWalkBurst, isWalking, cells, activeCell, walkPath
}) => {
  const [isCapturing, setIsCapturing] = useState(false);
  const [showKernelMenu, setShowKernelMenu] = useState(false);
  const stageLabels = ['Manifold', 'KNN Graph', 'Kernel', 'Matrix', 'Walk', 'Macrostates'];
  const selectedCell = activeCell !== null ? cells[activeCell] : null;
  const uniqueVisited = walkPath.length > 0 ? new Set(walkPath).size : 0;

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
        const probs = activeCell !== null ? getTransitionProbs(activeCell, cells, kernel) : [];
        
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
            <h2 className="text-2xl font-bold text-red-500">Terminal States</h2>
            <p className="text-slate-300">
              By analyzing the Markov Chain (using estimators like GPCCA+), CellRank automatically identifies <strong>Macrostates</strong>.
            </p>
            <p className="text-slate-300">
              These are the "sinks" of the process—states where random walks tend to end up. These correspond to the mature cell types.
            </p>
             <div className="mt-4 p-4 border border-red-500/30 bg-red-900/20 rounded">
                <h3 className="font-bold text-red-400 mb-2">Fate Probability</h3>
                <p className="text-sm">
                    <M>P</M>(fate<M><sub>A</sub></M> | cell<M><sub>i</sub></M>) is the probability that a random walk starting at <M>i</M> ends in Terminal State A.
                </p>
            </div>
           </div>
        )
      default:
        return null;
    }
  };

  return (
    <div className="absolute top-0 left-0 h-full w-full pointer-events-none flex flex-col justify-between">
      {/* Header */}
      <div className="bg-slate-900/90 backdrop-blur border-b border-slate-700 p-4 pointer-events-auto flex justify-between items-center">
        <div>
            <h1 className="text-xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-blue-400 to-emerald-400">
                CellRank 2 Explorer
            </h1>
            <p className="text-xs text-slate-400">Interactive Tutorial for Biologists</p>
        </div>
        <div className="flex items-center gap-4">
             <div className="hidden md:flex text-xs text-slate-500 gap-2">
                <span className="border border-slate-700 px-1 rounded bg-slate-800">←</span>
                <span className="border border-slate-700 px-1 rounded bg-slate-800">→</span>
                <span>to navigate</span>
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

      <div className="px-6 pt-3 pb-1 pointer-events-auto">
        <div className="max-w-3xl flex flex-wrap gap-2">
            {stageLabels.map((label, idx) => (
                <button
                    key={label}
                    onClick={() => setStage(idx as GameStage)}
                    className={`text-[11px] px-2.5 py-1 rounded-full border transition ${
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
      <div className="p-6 pointer-events-auto max-w-md">
        <div className="bg-slate-900/95 backdrop-blur border border-slate-700 p-6 rounded-xl shadow-2xl">
            {renderContent()}

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
      
      {/* Legend Footer */}
      <div className="p-4 bg-slate-900/80 backdrop-blur text-xs flex gap-6 justify-center text-slate-400 border-t border-slate-800">
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-yellow-400"></span> Stem</div>
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-slate-300"></span> Progenitor</div>
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-red-500"></span> Type A</div>
        <div className="flex items-center gap-2"><span className="w-3 h-3 rounded-full bg-blue-500"></span> Type B</div>
      </div>
    </div>
  );
};
