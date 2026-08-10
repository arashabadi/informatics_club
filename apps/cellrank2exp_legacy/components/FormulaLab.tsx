import React, { useState } from 'react';
import { ArrowLeft, Sigma, TrendingUp, Clock, Activity, ArrowRightLeft, GitMerge, Combine, FunctionSquare, Network } from 'lucide-react';

// ----------------------------------------------------------------------
// TYPES & NAVIGATION
// ----------------------------------------------------------------------

type ModuleId = 'COMPARE' | 'PSEUDOTIME' | 'VELOCITY' | 'CYTOTRACE' | 'REALTIME' | 'COMBINATION' | 'MATRIX' | 'GPCCA' | 'ABSORPTION';

const MODULES: { id: ModuleId; label: string; icon: React.ElementType }[] = [
    { id: 'COMPARE', label: 'Kernel Quick Compare', icon: FunctionSquare },
    { id: 'PSEUDOTIME', label: 'Pseudotime Kernel', icon: Clock },
    { id: 'VELOCITY', label: 'Velocity Kernel', icon: TrendingUp },
    { id: 'CYTOTRACE', label: 'CytoTRACE Kernel', icon: Activity },
    { id: 'REALTIME', label: 'RealTime (OT) Kernel', icon: ArrowRightLeft },
    { id: 'COMBINATION', label: 'Kernel Combination', icon: Combine },
    { id: 'MATRIX', label: 'Transition Matrix', icon: Sigma },
    { id: 'GPCCA', label: 'GPCCA & Macrostates', icon: GitMerge },
    { id: 'ABSORPTION', label: 'Absorption Probabilities', icon: Network },
];

export const FormulaLab: React.FC<{ onBack: () => void }> = ({ onBack }) => {
  const [activeModule, setActiveModule] = useState<ModuleId>('COMPARE');
  const activeIdx = MODULES.findIndex((m) => m.id === activeModule);
  const activeLabel = MODULES[activeIdx]?.label ?? '';

  return (
    <div className="w-full h-screen bg-slate-950 text-slate-100 font-sans flex overflow-hidden">
      {/* Sidebar */}
      <div className="w-64 bg-slate-900 border-r border-slate-800 flex flex-col">
        <div className="p-4 border-b border-slate-800 flex items-center gap-3">
             <button onClick={onBack} className="hover:bg-slate-800 p-2 rounded transition text-slate-400 hover:text-white">
                <ArrowLeft className="w-5 h-5" />
            </button>
            <h1 className="font-bold bg-clip-text text-transparent bg-gradient-to-r from-emerald-400 to-cyan-400">
                Formula Lab
            </h1>
        </div>
        <div className="flex-1 overflow-y-auto p-2 space-y-1">
            {MODULES.map((m) => (
                <button
                    key={m.id}
                    onClick={() => setActiveModule(m.id)}
                    className={`w-full px-4 py-3 rounded-lg text-sm font-bold flex items-center gap-3 transition-all ${
                        activeModule === m.id 
                        ? 'bg-blue-600/20 text-blue-400 border border-blue-500/50' 
                        : 'text-slate-400 hover:bg-slate-800 hover:text-slate-200'
                    }`}
                >
                    <m.icon size={16} /> {m.label}
                </button>
            ))}
        </div>
        <div className="p-4 border-t border-slate-800 text-xs text-slate-500">
            Hover over any math symbol to see its definition.
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 overflow-y-auto bg-slate-950 p-8">
        <div className="max-w-4xl mx-auto">
            <div className="mb-6 rounded-xl border border-slate-800 bg-slate-900/60 px-4 py-3 flex flex-wrap items-center justify-between gap-3">
                <div className="text-sm font-semibold text-slate-200">{activeLabel}</div>
                <div className="text-xs text-slate-400 font-mono">
                    Module {activeIdx + 1} / {MODULES.length}
                </div>
            </div>
            {activeModule === 'COMPARE' && <QuickCompareLab />}
            {activeModule === 'PSEUDOTIME' && <PseudotimeLab />}
            {activeModule === 'VELOCITY' && <VelocityLab />}
            {activeModule === 'CYTOTRACE' && <CytoTraceLab />}
            {activeModule === 'REALTIME' && <RealTimeLab />}
            {activeModule === 'COMBINATION' && <CombinationLab />}
            {activeModule === 'MATRIX' && <MatrixLab />}
            {activeModule === 'GPCCA' && <GpccaLab />}
            {activeModule === 'ABSORPTION' && <AbsorptionLab />}
        </div>
      </div>
    </div>
  );
};

// ----------------------------------------------------------------------
// SHARED COMPONENTS
// ----------------------------------------------------------------------

// Interactive Math Token
const T = ({ 
    val, 
    name, 
    desc 
}: { 
    val: React.ReactNode, 
    name: string, 
    desc: string 
}) => (
    <span className="group relative inline-block cursor-help border-b border-dotted border-slate-500 mx-1">
        <span className="font-serif italic font-bold text-slate-100">{val}</span>
        
        {/* Tooltip */}
        <span className="invisible group-hover:visible opacity-0 group-hover:opacity-100 transition-opacity absolute bottom-full left-1/2 -translate-x-1/2 mb-2 w-64 bg-slate-900 border border-slate-700 p-3 rounded-lg shadow-xl z-50 pointer-events-none">
            <span className="block text-emerald-400 font-bold text-xs mb-1 uppercase tracking-wider">{name}</span>
            <span className="block text-slate-300 text-xs leading-relaxed">{desc}</span>
            {/* Arrow */}
            <span className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-slate-900"></span>
        </span>
    </span>
);

const Section = ({ title, children }: { title: string, children?: React.ReactNode }) => (
    <div className="bg-slate-900/50 border border-slate-800 p-6 rounded-xl mb-6">
        <h2 className="text-lg font-bold text-slate-200 mb-6 flex items-center gap-2">
            {title}
        </h2>
        {children}
    </div>
);

const Range = ({ label, value, min, max, step, onChange, color = "blue" }: any) => (
    <div className="mb-4">
        <div className="flex justify-between text-sm mb-2 text-slate-400">
            <span>{label}</span>
            <span className="font-mono text-white">{value.toFixed(step < 1 ? 2 : 0)}</span>
        </div>
        <input 
            type="range" min={min} max={max} step={step}
            value={value} onChange={e => onChange(parseFloat(e.target.value))}
            className={`w-full h-2 rounded-lg appearance-none bg-slate-800 cursor-pointer accent-${color}-500`}
        />
    </div>
);

const QuickCompareLab = () => {
    const [pseudoDelta, setPseudoDelta] = useState(0.08);
    const [velocityCos, setVelocityCos] = useState(0.6);
    const [potencyDiff, setPotencyDiff] = useState(0.35);
    const [transportCost, setTransportCost] = useState(1.2);
    const [alpha, setAlpha] = useState(0.5);

    const pseudoRaw = Math.exp(8 * pseudoDelta);
    const velocityRaw = Math.exp(velocityCos / 0.35);
    const cytoRaw = 1 / (1 + Math.exp(-8 * potencyDiff));
    const realtimeRaw = Math.exp(-transportCost);
    const combinedRaw = alpha * velocityRaw + (1 - alpha) * pseudoRaw;

    const kernels = [
        {
            id: 'Pseudotime',
            raw: pseudoRaw,
            color: 'bg-blue-500',
            text: 'text-blue-300',
            hint: pseudoDelta >= 0 ? 'Forward pseudotime favored' : 'Backward edge penalized',
        },
        {
            id: 'Velocity',
            raw: velocityRaw,
            color: 'bg-purple-500',
            text: 'text-purple-300',
            hint: velocityCos >= 0 ? 'Aligned with velocity field' : 'Opposes velocity field',
        },
        {
            id: 'CytoTRACE',
            raw: cytoRaw,
            color: 'bg-yellow-500',
            text: 'text-yellow-300',
            hint: potencyDiff >= 0 ? 'High-to-low potency flow' : 'Potency increase discouraged',
        },
        {
            id: 'RealTime',
            raw: realtimeRaw,
            color: 'bg-cyan-500',
            text: 'text-cyan-300',
            hint: transportCost <= 2 ? 'Low OT transport cost' : 'High OT transport cost',
        },
        {
            id: 'Combined',
            raw: combinedRaw,
            color: 'bg-indigo-500',
            text: 'text-indigo-300',
            hint: `Blend α=${alpha.toFixed(2)} of Velocity + (1-α) Pseudotime`,
        },
    ];

    const total = kernels.reduce((acc, k) => acc + k.raw, 0);
    const normalized = kernels.map((k) => ({
        ...k,
        prob: total > 0 ? k.raw / total : 0,
    }));
    const winner = normalized.reduce((best, current) => (current.prob > best.prob ? current : best));

    return (
        <div className="space-y-6">
            <Section title="Kernel Quick Compare">
                <p className="text-sm text-slate-300 leading-relaxed">
                    Tune one shared scenario and compare how each kernel scores the same transition candidate.
                    This makes kernel assumptions easy to contrast before full transition-matrix construction.
                </p>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                <Section title="Shared Scenario Controls">
                    <Range
                        label="Pseudotime Difference (τ_j - τ_i)"
                        value={pseudoDelta}
                        min={-0.3}
                        max={0.3}
                        step={0.01}
                        onChange={setPseudoDelta}
                        color="blue"
                    />
                    <Range
                        label="Velocity Alignment cos(v_i, x_j-x_i)"
                        value={velocityCos}
                        min={-1}
                        max={1}
                        step={0.05}
                        onChange={setVelocityCos}
                        color="purple"
                    />
                    <Range
                        label="Potency Difference (pot_i - pot_j)"
                        value={potencyDiff}
                        min={-0.8}
                        max={0.8}
                        step={0.02}
                        onChange={setPotencyDiff}
                        color="yellow"
                    />
                    <Range
                        label="OT Transport Cost"
                        value={transportCost}
                        min={0}
                        max={6}
                        step={0.1}
                        onChange={setTransportCost}
                        color="cyan"
                    />
                    <Range
                        label="Combination Weight α (Velocity)"
                        value={alpha}
                        min={0}
                        max={1}
                        step={0.05}
                        onChange={setAlpha}
                        color="indigo"
                    />
                </Section>

                <Section title="Normalized Kernel Scores">
                    <div className="space-y-4">
                        {normalized.map((k) => (
                            <div key={k.id}>
                                <div className="flex justify-between items-center text-xs mb-1">
                                    <span className={`font-bold ${k.text}`}>{k.id}</span>
                                    <span className="font-mono text-slate-300">{(k.prob * 100).toFixed(1)}%</span>
                                </div>
                                <div className="h-3 rounded-full bg-slate-800 overflow-hidden">
                                    <div className={`h-full ${k.color}`} style={{ width: `${k.prob * 100}%` }} />
                                </div>
                                <div className="text-[11px] text-slate-500 mt-1">{k.hint}</div>
                            </div>
                        ))}
                    </div>
                    <div className="mt-5 p-3 rounded-lg border border-emerald-500/40 bg-emerald-900/20">
                        <div className="text-xs uppercase tracking-wider text-emerald-300 mb-1">Dominant Kernel</div>
                        <div className="text-sm font-semibold text-emerald-200">
                            {winner.id} drives this transition under the current settings.
                        </div>
                    </div>
                </Section>
            </div>
        </div>
    );
};

// ----------------------------------------------------------------------
// MODULE 1: Pseudotime Kernel
// ----------------------------------------------------------------------
const PseudotimeLab = () => {
    const [ti, setTi] = useState(0.5);
    const [tj, setTj] = useState(0.6);
    const [b, setB] = useState(20);

    const diff = tj - ti;
    const weight = Math.exp(b * diff);
    const displayWeight = weight > 100 ? '> 100' : weight.toFixed(4);

    return (
        <div className="space-y-6">
            <Section title="Formula Definition">
                <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <T val="P(j|i)" name="Transition Probability" desc="Probability of transitioning from cell i to cell j." />
                    <span>∝</span>
                    <span>exp(</span>
                    <T val="b" name="Kernel Strength" desc="A parameter that controls how strongly we favor transitions into the future. Higher values mean stricter forward movement." />
                    <span>· (</span>
                    <T val="τ_j" name="Pseudotime of Neighbor" desc="The estimated developmental progress (0 to 1) of the target neighbor cell." />
                    <span>-</span>
                    <T val="τ_i" name="Pseudotime of Cell" desc="The estimated developmental progress (0 to 1) of the current cell." />
                    <span>))</span>
                </div>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                <Section title="Parameters">
                    <Range label="Cell i Pseudotime (τ_i)" value={ti} min={0} max={1} step={0.01} onChange={setTi} color="blue" />
                    <Range label="Cell j Pseudotime (τ_j)" value={tj} min={0} max={1} step={0.01} onChange={setTj} color="pink" />
                    <Range label="Strength (b)" value={b} min={1} max={50} step={1} onChange={setB} color="emerald" />
                </Section>

                <Section title="Visualization">
                     <div className="h-40 relative flex items-center justify-center bg-slate-900 rounded-lg overflow-hidden">
                        <div className="absolute h-1 w-3/4 bg-slate-700 rounded-full" />
                        
                        {/* Cell i */}
                        <div className="absolute top-1/2 -translate-y-1/2 w-8 h-8 rounded-full bg-blue-500 border-2 border-white flex items-center justify-center z-10 transition-all" style={{ left: `${10 + ti * 70}%` }}>
                            <span className="font-serif italic text-xs font-bold">i</span>
                        </div>
                        
                        {/* Cell j */}
                        <div className="absolute top-1/2 -translate-y-1/2 w-8 h-8 rounded-full bg-pink-500 border-2 border-white flex items-center justify-center z-10 transition-all" style={{ left: `${10 + tj * 70}%` }}>
                            <span className="font-serif italic text-xs font-bold">j</span>
                        </div>

                         {/* Arrow */}
                         {Math.abs(diff) > 0.02 && (
                            <svg className="absolute inset-0 w-full h-full pointer-events-none">
                                <defs>
                                    <marker id="arrow" markerWidth="10" markerHeight="7" refX="0" refY="3.5" orient="auto">
                                        <polygon points="0 0, 10 3.5, 0 7" fill={diff > 0 ? "#4ade80" : "#ef4444"} />
                                    </marker>
                                </defs>
                                <line 
                                    x1={`${10 + ti * 70}%`} y1="50%" 
                                    x2={`${10 + tj * 70}%`} y2="50%" 
                                    stroke={diff > 0 ? "#4ade80" : "#ef4444"} 
                                    strokeWidth={Math.min(8, Math.max(1, weight / 5))} 
                                    markerEnd="url(#arrow)" 
                                />
                            </svg>
                         )}
                     </div>
                     <div className="mt-4 text-center font-mono text-xl">
                        Weight = <span className={diff > 0 ? "text-green-400" : "text-red-400"}>{displayWeight}</span>
                     </div>
                </Section>
            </div>
        </div>
    );
};

// ----------------------------------------------------------------------
// MODULE 2: Velocity Kernel
// ----------------------------------------------------------------------
const VelocityLab = () => {
    const [angle, setAngle] = useState(0);
    const [sigma, setSigma] = useState(0.5);

    // cos(angle)
    const cosSim = Math.cos(angle * Math.PI / 180);
    // Softmax-like weight: exp(cos / sigma)
    const weight = Math.exp(cosSim / sigma);

    return (
        <div className="space-y-6">
            <Section title="Formula Definition">
                <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <T val="P(j|i)" name="Transition Probability" desc="Probability of moving from cell i to cell j." />
                    <span>∝</span>
                    <span>exp(</span>
                    <T val="cos(v_i, x_j - x_i)" name="Cosine Similarity" desc="Measures alignment between the cell's internal velocity vector (v_i) and the direction towards the neighbor (x_j - x_i)." />
                    <span>/</span>
                    <T val="σ" name="Bandwidth (Sigma)" desc="Controls the 'width' of the kernel. Lower values make the kernel more 'peaky', only allowing transitions that perfectly align with velocity." />
                    <span>)</span>
                </div>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                <Section title="Parameters">
                    <Range label="Angle Difference (°)" value={angle} min={0} max={180} step={1} onChange={setAngle} color="purple" />
                    <Range label="Bandwidth (σ)" value={sigma} min={0.1} max={2.0} step={0.1} onChange={setSigma} color="emerald" />
                    <p className="text-xs text-slate-500 mt-2">
                        Angle 0° means the neighbor is exactly where the velocity vector points. 180° means it's behind.
                    </p>
                </Section>

                <Section title="Visualization">
                     <div className="h-48 relative flex items-center justify-center bg-slate-900 rounded-lg">
                        {/* Center i */}
                        <div className="w-4 h-4 bg-purple-500 rounded-full z-20 shadow-[0_0_10px_#a855f7]" />
                        
                        {/* Velocity Vector (Fixed Right) */}
                        <div className="absolute w-24 h-0.5 bg-purple-500 origin-left left-1/2 flex items-center">
                            <span className="absolute right-0 -top-4 text-purple-400 text-xs font-bold">v_i</span>
                            <div className="absolute right-0 w-2 h-2 border-t-2 border-r-2 border-purple-500 rotate-45 transform translate-x-1" />
                        </div>

                        {/* Neighbor Direction (Rotatable) */}
                        <div 
                            className="absolute w-24 h-0.5 bg-slate-600 origin-left left-1/2 border-t border-dashed border-slate-400"
                            style={{ transform: `rotate(${angle}deg)` }}
                        >
                            <div className="absolute right-0 w-4 h-4 bg-orange-500 rounded-full -translate-y-1/2 shadow shadow-orange-500/50 flex items-center justify-center">
                                <span className="text-[8px] font-bold text-black">j</span>
                            </div>
                        </div>

                        {/* Angle Arc */}
                        <svg className="absolute w-full h-full pointer-events-none">
                            <path 
                                d={`M ${50 + 10}% 50% A 40 40 0 0 1 ${50 + 10 * Math.cos(angle * Math.PI/180)}% ${50 + 10 * Math.sin(angle * Math.PI/180)}%`}
                                fill="none"
                                stroke="rgba(255,255,255,0.2)"
                                className={angle > 0 ? "visible" : "invisible"}
                            />
                        </svg>
                     </div>
                     <div className="mt-4 grid grid-cols-2 gap-4 text-center font-mono text-sm">
                        <div className="bg-slate-800 p-2 rounded">Cos Sim: <span className="text-white">{cosSim.toFixed(3)}</span></div>
                        <div className="bg-slate-800 p-2 rounded">Weight: <span className="text-emerald-400">{weight.toFixed(3)}</span></div>
                     </div>
                </Section>
            </div>
        </div>
    );
}

// ----------------------------------------------------------------------
// MODULE 3: CytoTRACE Kernel
// ----------------------------------------------------------------------
const CytoTraceLab = () => {
    const [genesI, setGenesI] = useState(5000); // Stem-like
    const [genesJ, setGenesJ] = useState(4000); // Diff-like

    // Simple logistic-like probability based on difference
    // CytoTRACE assumes: High Gene Count -> Low Gene Count
    const diff = genesI - genesJ; 
    // If I has more genes than J, flow is I->J (Prob high)
    // Formula approx: 1 / (1 + exp(-k * (score_i - score_j)))
    const prob = 1 / (1 + Math.exp(-0.005 * diff));

    return (
        <div className="space-y-6">
            <Section title="Formula Definition">
                <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <T val="P(j|i)" name="Transition Probability" desc="Probability flow based on developmental potential." />
                    <span>∝</span>
                    <T val="Potency(i)" name="Developmental Potency" desc="A score representing how 'stem-like' a cell is. In CytoTRACE, this is the Number of Expressed Genes." />
                    <span> {'>'} </span>
                    <T val="Potency(j)" name="Developmental Potency" desc="We expect differentiation to flow from High Potency to Low Potency." />
                </div>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                <Section title="Parameters">
                    <Range label="Expressed Genes (Cell i)" value={genesI} min={1000} max={10000} step={100} onChange={setGenesI} color="yellow" />
                    <Range label="Expressed Genes (Cell j)" value={genesJ} min={1000} max={10000} step={100} onChange={setGenesJ} color="blue" />
                    <p className="text-xs text-slate-500 mt-2">
                        Stem cells generally express more genes (transcriptional diversity). Differentiated cells silence unnecessary genes.
                    </p>
                </Section>

                <Section title="Visualization">
                     <div className="flex items-center justify-center gap-8 h-40 bg-slate-900 rounded-lg p-4">
                        {/* Cell I */}
                        <div className="flex flex-col items-center gap-2 relative">
                             <div 
                                className="rounded-full bg-yellow-500 border-2 border-white transition-all duration-500 flex items-center justify-center"
                                style={{ width: genesI / 100, height: genesI / 100, opacity: 0.8 }}
                             >
                                <span className="font-serif italic font-bold text-black">i</span>
                             </div>
                             <span className="text-xs text-yellow-500">{genesI} Genes</span>
                             <span className="text-[10px] text-slate-500 uppercase">Source</span>
                        </div>

                        {/* Arrow */}
                         <div className="flex-1 flex flex-col items-center">
                            <div className="h-1 bg-slate-700 w-full relative">
                                <div 
                                    className={`absolute top-0 left-0 h-full transition-all duration-300 ${prob > 0.5 ? 'bg-green-500' : 'bg-red-500'}`}
                                    style={{ width: `${prob * 100}%` }}
                                />
                            </div>
                            <span className={`mt-2 font-bold font-mono ${prob > 0.5 ? 'text-green-400' : 'text-red-400'}`}>
                                P(i→j) = {(prob * 100).toFixed(1)}%
                            </span>
                            <span className="text-xs text-slate-500">
                                {genesI > genesJ ? "Potency Loss (Expected)" : "Potency Gain (Rare)"}
                            </span>
                         </div>

                        {/* Cell J */}
                        <div className="flex flex-col items-center gap-2">
                             <div 
                                className="rounded-full bg-blue-500 border-2 border-white transition-all duration-500 flex items-center justify-center"
                                style={{ width: genesJ / 100, height: genesJ / 100, opacity: 0.8 }}
                             >
                                <span className="font-serif italic font-bold text-white">j</span>
                             </div>
                             <span className="text-xs text-blue-500">{genesJ} Genes</span>
                             <span className="text-[10px] text-slate-500 uppercase">Target</span>
                        </div>
                     </div>
                </Section>
            </div>
        </div>
    )
}

// ----------------------------------------------------------------------
// MODULE 4: RealTime (Optimal Transport)
// ----------------------------------------------------------------------
const RealTimeLab = () => {
    // 1D Optimal Transport visualization
    const [dist, setDist] = useState(2);
    const [epsilon, setEpsilon] = useState(0.5);

    // Cost = distance^2
    const cost = dist * dist;
    // Entropic regularization effect (simplified visual)
    // If epsilon is high, coupling is fuzzy (more spread). If low, strictly closest.
    const fuzziness = epsilon * 10; 

    return (
         <div className="space-y-6">
            <Section title="Formula Definition (Optimal Transport)">
                <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <T val="π*" name="Optimal Coupling" desc="The matrix describing how much mass moves from cell i at time t to cell j at time t+1." />
                    <span>= argmin</span>
                    <span>( Σ</span>
                    <T val="c(x_i, y_j)" name="Cost Function" desc="The cost to move mass. Usually squared Euclidean distance in gene expression space." />
                    <span>·</span>
                    <T val="π_ij" name="Transport Mass" desc="Amount of mass moved from i to j." />
                    <span>-</span>
                    <T val="ε" name="Epsilon" desc="Entropic regularization parameter. Higher values make the mapping 'fuzzier' and faster to compute." />
                    <T val="H(π)" name="Entropy" desc="Entropy of the coupling matrix. Encourages spreading mass to multiple similar neighbors." />
                    <span>)</span>
                </div>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                 <Section title="Parameters">
                    <Range label="Distance in Gene Space (Cost)" value={dist} min={0} max={10} step={0.5} onChange={setDist} color="red" />
                    <Range label="Regularization (ε)" value={epsilon} min={0.1} max={2.0} step={0.1} onChange={setEpsilon} color="cyan" />
                </Section>
                
                <Section title="Visualization">
                    <div className="h-40 bg-slate-900 rounded-lg p-4 relative flex items-center justify-center">
                        {/* Time T1 */}
                        <div className="absolute left-8 h-32 w-2 border-l border-dashed border-slate-600 flex flex-col items-center justify-center">
                            <span className="absolute -top-6 text-xs text-slate-500 w-20 text-center">Time t</span>
                            <div className="w-4 h-4 bg-yellow-500 rounded-full z-10"></div>
                        </div>

                         {/* Time T2 */}
                         <div className="absolute right-8 h-32 w-2 border-l border-dashed border-slate-600 flex flex-col items-center justify-center">
                            <span className="absolute -top-6 text-xs text-slate-500 w-20 text-center">Time t+1</span>
                            {/* Target Cell moves up/down based on 'distance' visual abstraction */}
                            <div className="w-4 h-4 bg-purple-500 rounded-full z-10" style={{ transform: `translateY(${dist * 10 - 20}px)` }}></div>
                        </div>

                        {/* Connection */}
                        <svg className="absolute inset-0 w-full h-full pointer-events-none">
                            {/* Main transport line */}
                            <line 
                                x1="40" y1="50%" 
                                x2="calc(100% - 40px)" y2={`${50 + (dist * 10 - 20) / 1.5}%`} 
                                stroke="#94a3b8" 
                                strokeWidth="2"
                            />
                            {/* Fuzziness aura */}
                            <line 
                                x1="40" y1="50%" 
                                x2="calc(100% - 40px)" y2={`${50 + (dist * 10 - 20) / 1.5}%`} 
                                stroke="rgba(56, 189, 248, 0.3)" 
                                strokeWidth={fuzziness}
                                strokeLinecap="round"
                            />
                        </svg>

                        <div className="absolute bottom-2 text-center w-full font-mono text-xs text-slate-400">
                            Transport Cost: {cost.toFixed(2)} | Fuzziness: {(fuzziness * 10).toFixed(0)}%
                        </div>
                    </div>
                </Section>
            </div>
         </div>
    )
}

// ----------------------------------------------------------------------
// MODULE 5: Kernel Combination
// ----------------------------------------------------------------------
const CombinationLab = () => {
    const [alpha, setAlpha] = useState(0.5);
    // Simulating two different kernels giving different transition probs for same edge
    const probK1 = 0.9; // e.g. Velocity says YES
    const probK2 = 0.2; // e.g. Similarity says NO

    const combined = alpha * probK1 + (1 - alpha) * probK2;

    return (
        <div className="space-y-6">
             <Section title="Formula Definition">
                <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <T val="T_comb" name="Combined Matrix" desc="The final transition matrix used for analysis." />
                    <span>=</span>
                    <T val="α" name="Weight" desc="User-defined weight parameter (0 to 1)." />
                    <span>·</span>
                    <T val="T_1" name="Kernel 1" desc="First kernel transition matrix (e.g. Velocity)." />
                    <span>+ (1 - α) ·</span>
                    <T val="T_2" name="Kernel 2" desc="Second kernel transition matrix (e.g. Connectivity)." />
                </div>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                 <Section title="Parameters">
                    <Range label="Combination Weight (α)" value={alpha} min={0} max={1} step={0.05} onChange={setAlpha} color="indigo" />
                    <div className="flex justify-between text-xs text-slate-500 px-1">
                        <span>100% Kernel 2</span>
                        <span>50/50</span>
                        <span>100% Kernel 1</span>
                    </div>
                </Section>

                <Section title="Visualization">
                    <div className="bg-slate-900 p-6 rounded-lg flex flex-col gap-4">
                        {/* Bar 1 */}
                        <div className="flex items-center gap-4">
                            <span className="w-20 text-xs font-bold text-pink-400">Kernel 1 (Velocity)</span>
                            <div className="flex-1 h-4 bg-slate-800 rounded-full overflow-hidden">
                                <div className="h-full bg-pink-500" style={{ width: `${probK1 * 100}%` }}></div>
                            </div>
                            <span className="w-8 text-xs">{probK1}</span>
                        </div>
                        
                        {/* Bar 2 */}
                        <div className="flex items-center gap-4">
                            <span className="w-20 text-xs font-bold text-blue-400">Kernel 2 (Similarity)</span>
                            <div className="flex-1 h-4 bg-slate-800 rounded-full overflow-hidden">
                                <div className="h-full bg-blue-500" style={{ width: `${probK2 * 100}%` }}></div>
                            </div>
                            <span className="w-8 text-xs">{probK2}</span>
                        </div>

                         {/* Result */}
                         <div className="mt-4 pt-4 border-t border-slate-700 flex items-center gap-4">
                            <span className="w-20 text-xs font-bold text-white">Combined</span>
                            <div className="flex-1 h-6 bg-slate-800 rounded-full overflow-hidden relative">
                                <div 
                                    className="h-full transition-all duration-300 bg-gradient-to-r from-blue-500 to-pink-500" 
                                    style={{ width: `${combined * 100}%` }}
                                ></div>
                            </div>
                            <span className="w-8 font-bold text-white">{combined.toFixed(2)}</span>
                        </div>
                    </div>
                </Section>
            </div>
        </div>
    )
}

// ----------------------------------------------------------------------
// MODULE 6: Transition Matrix (Existing but styled)
// ----------------------------------------------------------------------
const MatrixLab = () => {
    // 3 Neighbors with raw weights
    const [w1, setW1] = useState(5.0);
    const [w2, setW2] = useState(2.0);
    const [w3, setW3] = useState(1.0);

    const sum = w1 + w2 + w3;
    const p1 = w1 / sum;
    const p2 = w2 / sum;
    const p3 = w3 / sum;

    return (
        <div className="space-y-6">
             <Section title="Formula Definition">
                <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <T val="P_ij" name="Transition Probability" desc="Normalized probability of transitioning from i to j." />
                    <span>=</span>
                    <div className="flex flex-col items-center">
                        <span className="border-b border-slate-500 w-full text-center pb-1">
                            <T val="w_ij" name="Kernel Weight" desc="Raw similarity score from the kernel." />
                        </span>
                        <span className="pt-1">
                            <T val="Σ_k w_ik" name="Row Sum" desc="Sum of all outgoing weights from cell i. We divide by this to ensure probabilities sum to 1." />
                        </span>
                    </div>
                </div>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                 <Section title="Raw Weights (Kernel Output)">
                    <Range label="Neighbor 1 Weight" value={w1} min={0.1} max={10} step={0.1} onChange={setW1} color="red" />
                    <Range label="Neighbor 2 Weight" value={w2} min={0.1} max={10} step={0.1} onChange={setW2} color="yellow" />
                    <Range label="Neighbor 3 Weight" value={w3} min={0.1} max={10} step={0.1} onChange={setW3} color="blue" />
                </Section>

                <Section title="Row Stochastic Probabilities">
                     <div className="h-12 w-full bg-slate-800 rounded-lg flex overflow-hidden border border-slate-600">
                         <div style={{ width: `${p1 * 100}%` }} className="bg-red-500 h-full flex items-center justify-center text-white font-bold text-xs">{(p1*100).toFixed(0)}%</div>
                         <div style={{ width: `${p2 * 100}%` }} className="bg-yellow-500 h-full flex items-center justify-center text-black font-bold text-xs">{(p2*100).toFixed(0)}%</div>
                         <div style={{ width: `${p3 * 100}%` }} className="bg-blue-500 h-full flex items-center justify-center text-white font-bold text-xs">{(p3*100).toFixed(0)}%</div>
                    </div>
                    <p className="text-center mt-2 text-xs text-slate-500">Total Probability = 1.0</p>
                </Section>
            </div>
        </div>
    )
}

// ----------------------------------------------------------------------
// MODULE 7: GPCCA & Macrostates
// ----------------------------------------------------------------------
const GpccaLab = () => {
    const [gap, setGap] = useState(0.8); // Eigenvalue gap

    // Visualizing the Eigendecomposition
    // Top 3 eigenvalues: 1, 1, 0.99 (metastable) vs 1, 0.5, 0.4 (fast mixing)
    // If gap is high: 1, 1, gap...
    const eigs = [1.0, 0.98, gap, gap - 0.1, 0.1, 0.05];

    return (
        <div className="space-y-6">
            <Section title="Formula Definition (GPCCA)">
                 <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <span>X = </span>
                    <T val="arg max_X" name="Membership Matrix" desc="Matrix assigning cells to macrostates." />
                    <T val="Tr(X'TX)" name="Trace Optimization" desc="We look for a partition that maximizes the metastability (self-transition probability) of the macrostates." />
                    <span className="w-full text-center text-sm text-slate-500 mt-2">Based on Schur Decomposition: T = Q R Q'</span>
                </div>
            </Section>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
                <Section title="Spectrum Analysis">
                    <Range label="Spectral Gap" value={gap} min={0.2} max={0.9} step={0.05} onChange={setGap} color="orange" />
                    <p className="text-xs text-slate-400 mb-4">
                        Adjust the spectral gap. A large gap (low 3rd eigenvalue) indicates distinct, stable macrostates.
                    </p>
                    
                    {/* Eigenvalue Plot */}
                    <div className="h-40 bg-slate-900 border-l border-b border-slate-600 relative flex items-end justify-around px-4 pb-2">
                        {eigs.map((e, i) => (
                            <div key={i} className="flex flex-col items-center gap-1 group">
                                <div 
                                    className={`w-4 rounded-t ${i < 2 ? 'bg-emerald-400' : 'bg-slate-600'}`} 
                                    style={{ height: `${e * 100}%` }}
                                />
                                <span className="text-[10px] text-slate-500">λ_{i+1}</span>
                                <span className="absolute -top-6 text-xs bg-slate-800 p-1 rounded opacity-0 group-hover:opacity-100 transition-opacity">
                                    {e.toFixed(2)}
                                </span>
                            </div>
                        ))}
                    </div>
                </Section>

                <Section title="Macrostate Interpretation">
                     <div className="flex items-center justify-center h-48 gap-4">
                        {/* State A */}
                        <div className="w-20 h-20 rounded-full bg-red-500/20 border-2 border-red-500 flex items-center justify-center text-red-400 font-bold relative">
                            State A
                            {/* Self Loop */}
                             <svg className="absolute -top-4 -right-4 w-12 h-12 text-red-500">
                                <path d="M 0 10 C 10 -10, 30 -10, 20 10" fill="none" stroke="currentColor" strokeWidth="2" markerEnd="url(#arrow-red)" />
                            </svg>
                            <span className="absolute -top-6 -right-6 text-xs text-red-400">~0.99</span>
                        </div>
                        
                        {/* Mixing */}
                         <div className={`text-slate-500 font-mono transition-opacity ${gap > 0.7 ? 'opacity-20' : 'opacity-100'}`}>
                            &larr; Mixing &rarr;
                         </div>

                        {/* State B */}
                        <div className="w-20 h-20 rounded-full bg-blue-500/20 border-2 border-blue-500 flex items-center justify-center text-blue-400 font-bold relative">
                            State B
                            {/* Self Loop */}
                             <svg className="absolute -top-4 -right-4 w-12 h-12 text-blue-500">
                                <path d="M 0 10 C 10 -10, 30 -10, 20 10" fill="none" stroke="currentColor" strokeWidth="2" />
                            </svg>
                            <span className="absolute -top-6 -right-6 text-xs text-blue-400">~0.99</span>
                        </div>
                     </div>
                     <p className="text-center text-xs text-slate-400 mt-2">
                         {gap > 0.7 ? "High Metastability (Clear Terminal States)" : "Low Metastability (Transient/Cycling)"}
                     </p>
                </Section>
            </div>
        </div>
    )
}

// ----------------------------------------------------------------------
// MODULE 8: Absorption Probabilities
// ----------------------------------------------------------------------
const AbsorptionLab = () => {
    const [steps, setSteps] = useState(0); // For animation or logic

    // Simple 3 state Markov chain
    // 1 -> 2 (50%), 1 -> 1 (50%)
    // 2 -> 3 (Terminal)
    
    // Fate prob of Cell 1 reaching 3?
    // A_1 = P(1->2)*A_2 + P(1->1)*A_1
    // A_2 = 1.0
    // A_1 = 0.5 * 1 + 0.5 * A_1 => 0.5 A_1 = 0.5 => A_1 = 1.0 (eventually)
    
    return (
        <div className="space-y-6">
            <Section title="Formula Definition">
                 <div className="bg-slate-950 p-6 rounded-lg border border-slate-800 font-mono text-lg flex flex-wrap items-center justify-center gap-2">
                    <T val="a_i" name="Absorption Probability" desc="Probability that cell i will eventually end up in a specific terminal state." />
                    <span>=</span>
                    <T val="(I - T_{trans})^-1" name="Fundamental Matrix" desc="Inverse of (Identity - Transient Matrix). Represents expected number of visits to transient states." />
                    <span>·</span>
                    <T val="T_{trans->term}" name="Transition to Terminal" desc="Probabilities of moving from transient cells directly to terminal cells." />
                </div>
            </Section>

            <Section title="Example: Lineage Fate">
                <div className="relative h-40 bg-slate-900 rounded-lg p-8 flex items-center justify-between">
                    {/* Stem */}
                    <div className="flex flex-col items-center gap-2">
                        <div className="w-12 h-12 bg-yellow-500 rounded-full flex items-center justify-center text-black font-bold border-4 border-slate-800 z-10">
                            Stem
                        </div>
                         <div className="text-xs text-yellow-500">P(Fate A) = 50%</div>
                         <div className="text-xs text-yellow-500">P(Fate B) = 50%</div>
                    </div>

                    {/* Progenitors */}
                    <div className="absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2 w-full max-w-xs flex justify-between px-12">
                         <div className="w-8 h-8 bg-slate-600 rounded-full animate-pulse"></div>
                         <div className="w-8 h-8 bg-slate-600 rounded-full animate-pulse"></div>
                    </div>
                    
                    <svg className="absolute inset-0 w-full h-full pointer-events-none text-slate-600">
                        {/* Branches */}
                        <line x1="15%" y1="50%" x2="85%" y2="20%" stroke="currentColor" strokeWidth="2" strokeDasharray="4" />
                        <line x1="15%" y1="50%" x2="85%" y2="80%" stroke="currentColor" strokeWidth="2" strokeDasharray="4" />
                    </svg>

                    {/* Terminals */}
                    <div className="flex flex-col justify-between h-full">
                        <div className="flex items-center gap-2">
                            <div className="w-10 h-10 bg-red-500 rounded-full flex items-center justify-center text-white font-bold text-xs">A</div>
                            <span className="text-xs text-red-400">Fate A</span>
                        </div>
                        <div className="flex items-center gap-2">
                            <div className="w-10 h-10 bg-blue-500 rounded-full flex items-center justify-center text-white font-bold text-xs">B</div>
                            <span className="text-xs text-blue-400">Fate B</span>
                        </div>
                    </div>
                </div>
                <div className="mt-4 bg-slate-800 p-4 rounded text-sm text-slate-300">
                    The matrix inverse solves the linear system of equations to determine these probabilities globally, without needing to simulate millions of random walks.
                </div>
            </Section>
        </div>
    )
}
