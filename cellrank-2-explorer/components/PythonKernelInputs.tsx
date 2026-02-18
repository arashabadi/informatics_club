import React, { useMemo, useState } from 'react';
import { ArrowLeft, Check, Copy, Sparkles, Terminal } from 'lucide-react';

type SlotSnapshot = {
  slot: string;
  dtype: string;
  shape: string;
  example: string;
};

type ObsPreview = {
  columns: string[];
  rows: Array<Array<string | number>>;
};

type KernelInputCard = {
  id: string;
  label: string;
  emoji: string;
  colorClass: string;
  checklist: string[];
  funHint: string;
  structureSummary: string;
  slots: SlotSnapshot[];
  obsPreview: ObsPreview;
  matrixLabel: string;
  matrixPreview: number[][];
  snippet: string;
};

const formatPreviewValue = (value: string | number) => {
  if (typeof value === 'number') {
    return Number.isInteger(value) ? String(value) : value.toFixed(3);
  }
  return value;
};

const KERNEL_INPUTS: KernelInputCard[] = [
  {
    id: 'pseudotime',
    label: 'Pseudotime Kernel',
    emoji: '⏳',
    colorClass: 'from-emerald-500/20 to-blue-500/20 border-emerald-400/40',
    checklist: [
      '`adata.X` with cells × genes matrix',
      '`adata.obs["dpt_pseudotime"]` values in [0, 1]',
      'neighbors graph from `sc.pp.neighbors(adata)`',
    ],
    funHint: 'Challenge: perturb pseudotime by +0.05 and compare fate probabilities.',
    structureSummary: 'AnnData n_obs × n_vars = 3,698 × 2,000 (endocrine pancreas real dataset).',
    slots: [
      { slot: 'adata.X', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(3698, 2000)', example: 'log1p normalized expression counts' },
      { slot: 'adata.obs["dpt_pseudotime"]', dtype: 'pandas.Series<float64>', shape: '(3698,)', example: 'continuous lineage time from 0.0 to 1.0' },
      { slot: 'adata.obsp["connectivities"]', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(3698, 3698)', example: 'KNN graph (~15 non-zero neighbors/cell)' },
      { slot: 'adata.obsm["X_umap"]', dtype: 'numpy.ndarray<float32>', shape: '(3698, 2)', example: '2D embedding for plotting trajectories' },
    ],
    obsPreview: {
      columns: ['cell_id', 'cluster', 'dpt_pseudotime', 'batch'],
      rows: [
        ['cell_0001', 'Ngn3 low EP', 0.034, 'E14.5'],
        ['cell_0092', 'Ngn3 high EP', 0.281, 'E15.5'],
        ['cell_1403', 'Pre-endocrine', 0.572, 'E16.5'],
        ['cell_2669', 'Beta', 0.913, 'E16.5'],
      ],
    },
    matrixLabel: 'connectivities[0:4, 0:4]',
    matrixPreview: [
      [0.00, 0.12, 0.00, 0.03],
      [0.09, 0.00, 0.07, 0.00],
      [0.00, 0.11, 0.00, 0.16],
      [0.02, 0.00, 0.14, 0.00],
    ],
    snippet: `import scanpy as sc
import cellrank as cr

adata = sc.read_h5ad("pancreas.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.dpt(adata)

pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix()`,
  },
  {
    id: 'velocity',
    label: 'Velocity Kernel',
    emoji: '🧭',
    colorClass: 'from-fuchsia-500/20 to-rose-500/20 border-fuchsia-400/40',
    checklist: [
      '`adata.layers["spliced"]` and `adata.layers["unspliced"]`',
      '`adata.obsm["X_umap"]` or PCA for geometry',
      'RNA velocity estimated with `scvelo.tl.velocity`',
    ],
    funHint: 'Challenge: rerun with 30 vs 50 PCs and inspect lineage drift.',
    structureSummary: 'AnnData n_obs × n_vars = 2,531 × 2,000 with velocity layers and graph moments.',
    slots: [
      { slot: 'adata.layers["spliced"]', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(2531, 2000)', example: 'spliced counts per gene per cell' },
      { slot: 'adata.layers["unspliced"]', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(2531, 2000)', example: 'unspliced counts needed by dynamical model' },
      { slot: 'adata.obsm["velocity_umap"]', dtype: 'numpy.ndarray<float32>', shape: '(2531, 2)', example: 'velocity vectors projected onto UMAP' },
      { slot: 'adata.uns["velocity_graph"]', dtype: 'dict', shape: '{...}', example: 'stores graph metadata + parameters used by scVelo' },
    ],
    obsPreview: {
      columns: ['cell_id', 'cell_type', 'latent_time', 'velocity_confidence'],
      rows: [
        ['cell_0007', 'Ductal', 0.071, 0.821],
        ['cell_0840', 'Endocrine progenitor', 0.419, 0.744],
        ['cell_1301', 'Alpha', 0.778, 0.689],
        ['cell_1888', 'Beta', 0.904, 0.741],
      ],
    },
    matrixLabel: 'velocity_graph[0:4, 0:4]',
    matrixPreview: [
      [0.00, 0.21, 0.00, 0.05],
      [0.06, 0.00, 0.14, 0.00],
      [0.01, 0.09, 0.00, 0.18],
      [0.00, 0.02, 0.07, 0.00],
    ],
    snippet: `import scvelo as scv
import cellrank as cr

adata = scv.datasets.pancreas()
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()`,
  },
  {
    id: 'cytotrace',
    label: 'CytoTRACE Kernel',
    emoji: '🧬',
    colorClass: 'from-amber-500/20 to-orange-500/20 border-amber-400/40',
    checklist: [
      '`adata.X` with raw/normalized expression',
      '`adata.obs["ct_potency"]` where high = immature',
      'neighbors graph for local transitions',
    ],
    funHint: 'Challenge: compare top-driver genes for high vs low potency bins.',
    structureSummary: 'AnnData n_obs × n_vars = 5,212 × 2,500 with a potency score for each cell.',
    slots: [
      { slot: 'adata.X', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(5212, 2500)', example: 'expression matrix used for gene-count complexity' },
      { slot: 'adata.obs["ct_potency"]', dtype: 'pandas.Series<float64>', shape: '(5212,)', example: 'higher = more stem-like/transcriptionally diverse' },
      { slot: 'adata.var["highly_variable"]', dtype: 'pandas.Series<bool>', shape: '(2500,)', example: 'feature mask for stable graph construction' },
      { slot: 'adata.obsp["distances"]', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(5212, 5212)', example: 'pairwise neighborhood distances' },
    ],
    obsPreview: {
      columns: ['cell_id', 'lineage', 'ct_potency', 'n_genes_by_counts'],
      rows: [
        ['cell_0012', 'HSC', 0.944, 3891],
        ['cell_0264', 'MPP', 0.728, 3123],
        ['cell_1185', 'Myeloid prog.', 0.411, 2455],
        ['cell_4430', 'Neutrophil', 0.107, 1299],
      ],
    },
    matrixLabel: 'distance_graph[0:4, 0:4]',
    matrixPreview: [
      [0.00, 0.34, 0.00, 0.00],
      [0.31, 0.00, 0.18, 0.04],
      [0.00, 0.23, 0.00, 0.22],
      [0.00, 0.06, 0.19, 0.00],
    ],
    snippet: `import scanpy as sc
import cellrank as cr

adata = sc.read_h5ad("hematopoiesis.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.neighbors(adata)

# Example potency proxy: genes detected per cell
adata.obs["ct_potency"] = (adata.X > 0).sum(1).A1

ck = cr.kernels.CytoTRACEKernel(adata, potency_key="ct_potency")
ck.compute_transition_matrix()`,
  },
  {
    id: 'realtime',
    label: 'RealTime Kernel (OT)',
    emoji: '🕒',
    colorClass: 'from-cyan-500/20 to-sky-500/20 border-cyan-400/40',
    checklist: [
      '`adata.obs["time"]` with real sampling times',
      'cells from multiple timepoints (e.g., day0/day2/day5)',
      'consistent feature space across timepoints',
    ],
    funHint: 'Challenge: remove one intermediate timepoint and see transport changes.',
    structureSummary: 'AnnData n_obs × n_vars = 4,104 × 3,000 sampled across real experiment timepoints.',
    slots: [
      { slot: 'adata.obs["time"]', dtype: 'pandas.Categorical', shape: '(4104,)', example: 'ordered categories: day0 < day2 < day5 < day7' },
      { slot: 'adata.obsm["X_pca"]', dtype: 'numpy.ndarray<float32>', shape: '(4104, 50)', example: 'shared latent space used for OT cost' },
      { slot: 'adata.uns["time_order"]', dtype: 'list[str]', shape: '(4,)', example: 'explicit temporal ordering used by RealTimeKernel' },
      { slot: 'adata.obsp["connectivities"]', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(4104, 4104)', example: 'intra-timepoint neighborhood graph' },
    ],
    obsPreview: {
      columns: ['cell_id', 'time', 'cell_state', 'replicate'],
      rows: [
        ['cell_0021', 'day0', 'Neural stem', 'R1'],
        ['cell_0605', 'day2', 'Transit amplifying', 'R1'],
        ['cell_2280', 'day5', 'Neuroblast', 'R2'],
        ['cell_3794', 'day7', 'Neuron', 'R2'],
      ],
    },
    matrixLabel: 'OT coupling (local block)',
    matrixPreview: [
      [0.00, 0.17, 0.00, 0.00],
      [0.10, 0.00, 0.09, 0.00],
      [0.00, 0.12, 0.00, 0.15],
      [0.00, 0.00, 0.11, 0.00],
    ],
    snippet: `import scanpy as sc
import cellrank as cr

adata = sc.read_h5ad("timecourse_brain.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)

rk = cr.kernels.RealTimeKernel(adata, time_key="time")
rk.compute_transition_matrix()`,
  },
  {
    id: 'combined',
    label: 'Combined Kernel',
    emoji: '🎛️',
    colorClass: 'from-violet-500/20 to-indigo-500/20 border-violet-400/40',
    checklist: [
      'precomputed `VelocityKernel` + `PseudotimeKernel`',
      'matching `adata` reference for both kernels',
      'blend weight `alpha` for biological trust balance',
    ],
    funHint: 'Challenge: sweep α from 0.2 to 0.8 and plot fate stability.',
    structureSummary: 'Same AnnData object feeds two kernels; combined operator keeps a shared cell index space.',
    slots: [
      { slot: 'adata.obs["dpt_pseudotime"]', dtype: 'pandas.Series<float64>', shape: '(2531,)', example: 'global progress ordering from diffusion pseudotime' },
      { slot: 'adata.layers["spliced"/"unspliced"]', dtype: 'scipy.sparse.csr_matrix<float32>', shape: '(2531, 2000)', example: 'velocity features from RNA splicing kinetics' },
      { slot: 'vk.transition_matrix', dtype: 'scipy.sparse.csr_matrix<float64>', shape: '(2531, 2531)', example: 'velocity transition operator' },
      { slot: 'pk.transition_matrix', dtype: 'scipy.sparse.csr_matrix<float64>', shape: '(2531, 2531)', example: 'pseudotime transition operator' },
    ],
    obsPreview: {
      columns: ['cell_id', 'dpt_pseudotime', 'latent_time', 'cell_type'],
      rows: [
        ['cell_0044', 0.082, 0.091, 'Ductal'],
        ['cell_0914', 0.326, 0.447, 'Endocrine prog.'],
        ['cell_1607', 0.669, 0.701, 'Alpha'],
        ['cell_2210', 0.902, 0.878, 'Beta'],
      ],
    },
    matrixLabel: 'combined T = α·Tv + (1-α)·Tp',
    matrixPreview: [
      [0.00, 0.19, 0.00, 0.02],
      [0.07, 0.00, 0.10, 0.00],
      [0.00, 0.08, 0.00, 0.20],
      [0.01, 0.00, 0.11, 0.00],
    ],
    snippet: `import cellrank as cr

vk = cr.kernels.VelocityKernel(adata).compute_transition_matrix()
pk = cr.kernels.PseudotimeKernel(
    adata, time_key="dpt_pseudotime"
).compute_transition_matrix()

alpha = 0.6
combo = alpha * vk + (1 - alpha) * pk
combo.compute_transition_matrix()`,
  },
];

export const PythonKernelInputs: React.FC<{ onBack: () => void }> = ({ onBack }) => {
  const [selectedId, setSelectedId] = useState(KERNEL_INPUTS[0].id);
  const [copiedId, setCopiedId] = useState<string | null>(null);

  const selected = useMemo(
    () => KERNEL_INPUTS.find((item) => item.id === selectedId) ?? KERNEL_INPUTS[0],
    [selectedId]
  );

  const handleCopy = async () => {
    try {
      await navigator.clipboard.writeText(selected.snippet);
      setCopiedId(selected.id);
      window.setTimeout(() => setCopiedId((prev) => (prev === selected.id ? null : prev)), 1400);
    } catch {
      setCopiedId(null);
    }
  };

  return (
    <div className="w-full min-h-screen bg-slate-950 text-white p-4 md:p-6">
      <div className="max-w-7xl mx-auto">
        <div className="flex items-center justify-between gap-3 mb-4">
          <button
            onClick={onBack}
            className="inline-flex items-center gap-2 rounded-lg border border-slate-700 bg-slate-900 px-3 py-2 text-sm hover:bg-slate-800"
          >
            <ArrowLeft className="w-4 h-4" />
            Back to Menu
          </button>
          <div className="inline-flex items-center gap-2 text-xs text-slate-300 border border-slate-700 bg-slate-900/90 rounded-lg px-3 py-2">
            <Sparkles className="w-4 h-4 text-emerald-400" />
            Real data inputs + runnable CellRank kernel setup
          </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-[280px,1fr] gap-4">
          <div className="rounded-xl border border-slate-800 bg-slate-900 p-3 space-y-2 h-fit lg:sticky lg:top-4">
            {KERNEL_INPUTS.map((kernel) => (
              <button
                key={kernel.id}
                onClick={() => setSelectedId(kernel.id)}
                className={`w-full text-left rounded-lg border px-3 py-2 transition ${
                  selected.id === kernel.id
                    ? 'border-blue-400/50 bg-blue-500/10'
                    : 'border-slate-700 bg-slate-950 hover:bg-slate-800/70'
                }`}
              >
                <div className="text-sm font-semibold">{kernel.emoji} {kernel.label}</div>
              </button>
            ))}
          </div>

          <div className={`rounded-xl border bg-gradient-to-br ${selected.colorClass} p-4 md:p-6`}>
            <div className="flex flex-wrap items-center justify-between gap-3 mb-4">
              <h2 className="text-2xl font-bold">{selected.emoji} {selected.label}</h2>
              <button
                onClick={handleCopy}
                className="inline-flex items-center gap-2 rounded-lg border border-slate-600 bg-slate-900/80 px-3 py-2 text-sm hover:bg-slate-800"
              >
                {copiedId === selected.id ? <Check className="w-4 h-4 text-emerald-400" /> : <Copy className="w-4 h-4" />}
                {copiedId === selected.id ? 'Copied' : 'Copy Python'}
              </button>
            </div>

            <div className="grid grid-cols-1 xl:grid-cols-3 gap-4">
              <div className="rounded-lg border border-slate-700 bg-slate-900/90 p-4">
                <div className="flex items-center gap-2 text-sm text-emerald-300 mb-2">
                  <Terminal className="w-4 h-4" />
                  Input checklist
                </div>
                <ul className="space-y-2 text-sm text-slate-200">
                  {selected.checklist.map((line) => (
                    <li key={line}>• {line}</li>
                  ))}
                </ul>
                <div className="mt-4 rounded border border-emerald-500/30 bg-emerald-950/30 px-3 py-2 text-sm text-emerald-200">
                  {selected.funHint}
                </div>
              </div>

              <div className="rounded-lg border border-slate-700 bg-slate-900/90 p-4">
                <div className="text-xs uppercase tracking-wider text-slate-300 mb-2">Actual data structure snapshot</div>
                <p className="text-xs text-slate-300 mb-3">{selected.structureSummary}</p>

                <div className="space-y-2">
                  {selected.slots.map((slot) => (
                    <div key={slot.slot} className="rounded border border-slate-700 bg-slate-950/70 px-3 py-2">
                      <div className="flex items-center justify-between gap-2">
                        <span className="text-[11px] text-cyan-300 font-mono">{slot.slot}</span>
                        <span className="text-[11px] text-slate-500 font-mono">{slot.shape}</span>
                      </div>
                      <div className="text-[11px] text-slate-300 mt-0.5">{slot.dtype}</div>
                      <div className="text-[11px] text-slate-500">{slot.example}</div>
                    </div>
                  ))}
                </div>

                <div className="mt-3 overflow-x-auto">
                  <div className="text-[11px] uppercase tracking-wider text-slate-400 mb-1">adata.obs preview</div>
                  <table className="min-w-full text-[11px] border border-slate-700">
                    <thead className="bg-slate-900">
                      <tr>
                        {selected.obsPreview.columns.map((column) => (
                          <th key={column} className="px-2 py-1 text-left text-slate-300 border-b border-slate-700 whitespace-nowrap">{column}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {selected.obsPreview.rows.map((row, rowIndex) => (
                        <tr key={rowIndex} className="odd:bg-slate-900/30">
                          {row.map((cellValue, cellIndex) => (
                            <td key={cellIndex} className="px-2 py-1 text-slate-200 border-b border-slate-800 whitespace-nowrap">
                              {formatPreviewValue(cellValue)}
                            </td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>

                <div className="mt-3">
                  <div className="text-[11px] uppercase tracking-wider text-slate-400 mb-1">{selected.matrixLabel}</div>
                  <div
                    className="inline-grid gap-[2px] rounded border border-slate-700 bg-slate-950 p-2"
                    style={{ gridTemplateColumns: `repeat(${selected.matrixPreview[0]?.length ?? 0}, minmax(16px, 1fr))` }}
                  >
                    {selected.matrixPreview.map((row, rowIndex) =>
                      row.map((value, colIndex) => (
                        <div
                          key={`${rowIndex}-${colIndex}`}
                          className="w-6 h-6 rounded-[2px] text-[9px] flex items-center justify-center text-slate-100"
                          style={{ backgroundColor: `rgba(56, 189, 248, ${0.1 + value})` }}
                          title={`${value.toFixed(3)}`}
                        >
                          {value.toFixed(2)}
                        </div>
                      ))
                    )}
                  </div>
                </div>
              </div>

              <div className="rounded-lg border border-slate-700 bg-slate-950/95 p-4">
                <div className="text-xs uppercase tracking-wider text-slate-400 mb-2">Python recipe</div>
                <pre className="text-xs md:text-sm leading-relaxed overflow-auto text-slate-100 font-mono">
{selected.snippet}
                </pre>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};
