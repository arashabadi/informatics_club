import React from 'react';
import {
  ArrowLeft,
  ArrowRight,
  CheckCircle2,
  CircleDashed,
  FlaskConical,
  Sigma,
  GitMerge,
  Target,
  Activity,
} from 'lucide-react';

type FlowStep = {
  label: string;
  optional?: boolean;
  detail: string;
};

type FlowSection = {
  letter: 'a' | 'b' | 'c';
  title: string;
  icon: React.ElementType;
  required: FlowStep[];
  optional: FlowStep[];
};

const SECTION_STYLES: Record<FlowSection['letter'], { bar: string; chip: string; required: string }> = {
  a: {
    bar: 'from-sky-500/80 to-cyan-500/80',
    chip: 'text-sky-300',
    required: 'from-sky-900/70 to-cyan-900/60 border-sky-400/50',
  },
  b: {
    bar: 'from-indigo-500/80 to-fuchsia-500/80',
    chip: 'text-indigo-300',
    required: 'from-indigo-900/70 to-fuchsia-900/60 border-indigo-400/50',
  },
  c: {
    bar: 'from-pink-500/80 to-rose-500/80',
    chip: 'text-pink-300',
    required: 'from-pink-900/70 to-rose-900/60 border-pink-400/50',
  },
};

const SECTIONS: FlowSection[] = [
  {
    letter: 'a',
    title: 'Kernel',
    icon: FlaskConical,
    required: [
      { label: 'Initialization', detail: 'Load data view + kernel inputs.' },
      { label: 'Transition matrix', detail: 'Build row-stochastic transition matrix T.' },
    ],
    optional: [
      { label: 'Kernel combination', optional: true, detail: 'Fuse multiple kernels into one T.' },
      { label: 'Random walk simulation', optional: true, detail: 'Simulate trajectories over T.' },
      { label: 'CBC score', optional: true, detail: 'Cross-boundary correctness benchmark.' },
    ],
  },
  {
    letter: 'b',
    title: 'Estimator',
    icon: Sigma,
    required: [
      { label: 'Schur decomposition', detail: 'Decompose transition dynamics.' },
      { label: 'Macrostate computation', detail: 'Coarse-grain states with GPCCA.' },
      { label: 'Terminal state identification', detail: 'Detect absorbing/terminal states.' },
      { label: 'Fate probabilities', detail: 'Compute absorption probabilities.' },
    ],
    optional: [{ label: 'TSI score', optional: true, detail: 'Terminal-state identification score.' }],
  },
  {
    letter: 'c',
    title: 'Analysis',
    icon: Activity,
    required: [],
    optional: [
      { label: 'Lineage correlation', optional: true, detail: 'Correlate genes with lineage fate.' },
      { label: 'Putative driver ranking', optional: true, detail: 'Rank candidate lineage drivers.' },
      { label: 'GEX trends', optional: true, detail: 'Model lineage-specific expression trends.' },
      { label: 'Kernel comparison', optional: true, detail: 'Compare kernels across metrics.' },
    ],
  },
];

const StepCard: React.FC<{
  step: FlowStep;
  requiredTone: string;
}> = ({ step, requiredTone }) => {
  const cardClass = step.optional
    ? 'border-dashed border-slate-500/70 bg-slate-900/70'
    : `bg-gradient-to-r ${requiredTone}`;

  return (
    <div
      className={`min-h-[88px] rounded-xl border-2 p-4 flex flex-col justify-center gap-1 shadow-lg ${cardClass}`}
      title={step.optional ? 'Optional step' : 'Required step'}
    >
      <div className="flex items-center gap-2">
        {step.optional ? (
          <CircleDashed className="w-4 h-4 text-slate-300" />
        ) : (
          <CheckCircle2 className="w-4 h-4 text-emerald-300" />
        )}
        <span className="font-semibold text-slate-100">{step.label}</span>
      </div>
      <p className="text-xs text-slate-300">{step.detail}</p>
    </div>
  );
};

export const ProtocolWorkflow: React.FC<{ onBack: () => void }> = ({ onBack }) => {
  return (
    <div className="w-full min-h-screen bg-slate-950 text-slate-100">
      <div className="sticky top-0 z-20 bg-slate-900/95 border-b border-slate-800 backdrop-blur">
        <div className="max-w-7xl mx-auto px-6 py-4 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <button
              onClick={onBack}
              className="hover:bg-slate-800 p-2 rounded transition text-slate-400 hover:text-white"
              aria-label="Back to menu"
            >
              <ArrowLeft className="w-5 h-5" />
            </button>
            <div>
              <h1 className="text-xl font-bold text-white">CellRank Workflow Map</h1>
              <p className="text-xs text-slate-400">
                Full protocol view with required and optional steps
              </p>
            </div>
          </div>
          <div className="hidden sm:flex items-center gap-4 text-xs">
            <div className="flex items-center gap-2 text-slate-300">
              <span className="w-3 h-3 rounded-full bg-emerald-400" />
              Required step
            </div>
            <div className="flex items-center gap-2 text-slate-300">
              <span className="w-3 h-3 rounded-full border border-dashed border-slate-300" />
              Optional step
            </div>
          </div>
        </div>
      </div>

      <div className="max-w-7xl mx-auto px-6 py-8 space-y-6">
        {SECTIONS.map((section) => {
          const styles = SECTION_STYLES[section.letter];
          const Icon = section.icon;

          return (
            <section
              key={section.letter}
              className="rounded-2xl border border-slate-800 bg-slate-900/50 p-5 md:p-6"
            >
              <div className="flex items-center gap-4 mb-5">
                <div className={`w-3 h-14 rounded-full bg-gradient-to-b ${styles.bar}`} />
                <div className="flex items-center gap-3">
                  <span className={`text-4xl font-black leading-none ${styles.chip}`}>{section.letter}</span>
                  <div className="flex items-center gap-2">
                    <Icon className={`w-5 h-5 ${styles.chip}`} />
                    <h2 className="text-2xl font-bold text-white">{section.title}</h2>
                  </div>
                </div>
              </div>

              {section.required.length > 0 && (
                <div className="grid grid-cols-1 lg:grid-cols-[1fr_auto_1fr_auto_1fr_auto_1fr] gap-3 items-stretch mb-4">
                  {section.required.map((step, index) => (
                    <React.Fragment key={step.label}>
                      <StepCard step={step} requiredTone={styles.required} />
                      {index < section.required.length - 1 && (
                        <div className="hidden lg:flex items-center justify-center text-slate-500">
                          <ArrowRight className="w-5 h-5" />
                        </div>
                      )}
                    </React.Fragment>
                  ))}
                </div>
              )}

              {section.optional.length > 0 && (
                <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-4 gap-3">
                  {section.optional.map((step) => (
                    <StepCard key={step.label} step={step} requiredTone={styles.required} />
                  ))}
                </div>
              )}
            </section>
          );
        })}
      </div>
    </div>
  );
};
