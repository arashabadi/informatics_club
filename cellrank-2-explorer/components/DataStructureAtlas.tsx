import React, { useMemo, useState } from 'react';
import { ArrowLeft, ArrowRightLeft, Database, Sparkles } from 'lucide-react';

type SlotId = 'x' | 'layers' | 'obs' | 'var' | 'obsm' | 'varm' | 'obsp' | 'varp' | 'uns' | 'raw';

type SlotMap = {
  id: SlotId;
  annLabel: string;
  annPath: string;
  annMeaning: string;
  seuLabel: string;
  seuPath: string;
  seuMeaning: string;
  shapeHint: string;
  tip: string;
  pythonPeek: string;
  rPeek: string;
};

const SLOT_MAP: SlotMap[] = [
  {
    id: 'x',
    annLabel: 'X',
    annPath: 'adata.X',
    annMeaning: 'Primary expression matrix used for many downstream operations.',
    seuLabel: 'RNA data/counts',
    seuPath: 'obj[["RNA"]]$data / $counts',
    seuMeaning: 'Stored as assay layers in Seurat v5.',
    shapeHint: 'AnnData: cells × genes; Seurat assay layer: genes × cells',
    tip: 'Orientation is the biggest mental shift between AnnData and Seurat matrices.',
    pythonPeek: 'adata.X[:3, :5]',
    rPeek: 'LayerData(obj, assay = "RNA", layer = "data")[1:5, 1:3]',
  },
  {
    id: 'layers',
    annLabel: 'layers',
    annPath: 'adata.layers["counts"/"log1p"]',
    annMeaning: 'Alternative matrices aligned with X.',
    seuLabel: 'Layers()',
    seuPath: 'Layers(obj[["RNA"]])',
    seuMeaning: 'Assay5 is layer-native; this is the closest match to AnnData layers.',
    shapeHint: 'Same dimensions as X for each layer',
    tip: 'AnnData layers ≈ Seurat Assay5 layers.',
    pythonPeek: 'adata.layers["counts"][:3, :5]',
    rPeek: 'LayerData(obj, assay = "RNA", layer = "counts")[1:5, 1:3]',
  },
  {
    id: 'obs',
    annLabel: 'obs',
    annPath: 'adata.obs',
    annMeaning: 'Per-cell metadata table.',
    seuLabel: 'meta.data',
    seuPath: 'obj[[]] / obj@meta.data',
    seuMeaning: 'Cell-level metadata in Seurat object.',
    shapeHint: 'Rows = cells',
    tip: 'obs = observations = cells.',
    pythonPeek: 'adata.obs.head()',
    rPeek: 'head(obj[[]])',
  },
  {
    id: 'var',
    annLabel: 'var',
    annPath: 'adata.var',
    annMeaning: 'Per-gene metadata table.',
    seuLabel: 'assay feature meta',
    seuPath: 'obj[["RNA"]]@meta.data',
    seuMeaning: 'Feature-level metadata lives in assay metadata.',
    shapeHint: 'Rows = genes/features',
    tip: 'var = variables = genes/features.',
    pythonPeek: 'adata.var.head()',
    rPeek: 'head(obj[["RNA"]]@meta.data)',
  },
  {
    id: 'obsm',
    annLabel: 'obsm',
    annPath: 'adata.obsm["X_umap"/"X_pca"]',
    annMeaning: 'Per-cell multidimensional arrays (embeddings).',
    seuLabel: 'reductions',
    seuPath: 'Embeddings(obj[["umap"]])',
    seuMeaning: 'Stored in DimReduc objects.',
    shapeHint: 'Rows = cells, columns = embedding dimensions',
    tip: 'obsm is where most cell embeddings live.',
    pythonPeek: 'adata.obsm["X_umap"][:5, :2]',
    rPeek: 'Embeddings(obj[["umap"]])[1:5, 1:2]',
  },
  {
    id: 'varm',
    annLabel: 'varm',
    annPath: 'adata.varm["PCs"]',
    annMeaning: 'Per-gene multidimensional arrays (loadings).',
    seuLabel: 'feature loadings',
    seuPath: 'Loadings(obj[["pca"]])',
    seuMeaning: 'PCA loadings in DimReduc feature.loadings.',
    shapeHint: 'Rows = genes, columns = components',
    tip: 'varm mirrors obsm, but for genes.',
    pythonPeek: 'adata.varm["PCs"][:5, :3]',
    rPeek: 'Loadings(obj[["pca"]])[1:5, 1:3]',
  },
  {
    id: 'obsp',
    annLabel: 'obsp',
    annPath: 'adata.obsp["connectivities"/"distances"]',
    annMeaning: 'Cell-cell pairwise sparse matrices.',
    seuLabel: 'graphs/neighbors',
    seuPath: 'obj@graphs$RNA_snn / obj@neighbors',
    seuMeaning: 'KNN/SNN graph structures.',
    shapeHint: 'Square matrix over cells',
    tip: 'obsp is critical for graph-based trajectory methods.',
    pythonPeek: 'adata.obsp["connectivities"][:5, :5].A',
    rPeek: 'obj@graphs$RNA_snn[1:5, 1:5]',
  },
  {
    id: 'varp',
    annLabel: 'varp',
    annPath: 'adata.varp[...]',
    annMeaning: 'Gene-gene pairwise matrices.',
    seuLabel: 'custom misc/tools',
    seuPath: 'obj[["RNA"]]@misc$gene_graph',
    seuMeaning: 'No strict dedicated container; often project-defined.',
    shapeHint: 'Square matrix over genes',
    tip: 'varp exists explicitly in AnnData; Seurat usually needs conventions.',
    pythonPeek: 'adata.varp["gene_corr"][:5, :5]',
    rPeek: 'obj[["RNA"]]@misc$gene_graph[1:5, 1:5]',
  },
  {
    id: 'uns',
    annLabel: 'uns',
    annPath: 'adata.uns',
    annMeaning: 'Unstructured metadata/parameters.',
    seuLabel: 'misc/tools/commands',
    seuPath: 'obj@misc, obj@tools, obj@commands',
    seuMeaning: 'General container for settings, outputs, and run history.',
    shapeHint: 'Key-value tree',
    tip: 'uns is a metadata catch-all.',
    pythonPeek: 'list(adata.uns.keys())[:8]',
    rPeek: 'names(obj@misc)',
  },
  {
    id: 'raw',
    annLabel: 'raw',
    annPath: 'adata.raw.X + adata.raw.var',
    annMeaning: 'Frozen matrix + feature metadata snapshot.',
    seuLabel: 'raw counts layer/assay',
    seuPath: 'obj[["RNA"]]$counts or obj[["RNA_raw"]]',
    seuMeaning: 'Usually stored by assay/layer convention.',
    shapeHint: 'Snapshot, typically pre-transformation',
    tip: 'AnnData raw is explicit; Seurat raw is convention-based.',
    pythonPeek: 'adata.raw[:, :5].X[:3, :5]',
    rPeek: 'obj[["RNA"]]$counts[1:5, 1:3]',
  },
];

const ANN_LAYOUT: Record<SlotId, string> = {
  raw: 'col-start-1 row-start-1',
  var: 'col-start-3 col-span-2 row-start-1',
  obsp: 'col-start-1 row-start-2 row-span-2',
  obsm: 'col-start-2 row-start-2 row-span-2',
  x: 'col-start-3 row-start-2 row-span-2',
  obs: 'col-start-4 row-start-2 row-span-2',
  layers: 'col-start-2 col-span-2 row-start-4',
  varm: 'col-start-3 row-start-5',
  varp: 'col-start-4 row-start-5',
  uns: 'col-start-1 col-span-2 row-start-5',
};

const SEU_LAYOUT: Record<SlotId, string> = {
  raw: 'col-start-3 row-start-1',
  var: 'col-start-1 row-start-1',
  layers: 'col-start-2 row-start-1',
  obs: 'col-start-1 row-start-2',
  x: 'col-start-2 row-start-2',
  obsm: 'col-start-3 row-start-2',
  varm: 'col-start-1 row-start-3',
  obsp: 'col-start-2 row-start-3',
  varp: 'col-start-3 row-start-3',
  uns: 'col-start-1 col-span-3 row-start-4',
};

const colorFor = (id: SlotId) => {
  if (id === 'x' || id === 'layers') return 'from-emerald-500/30 to-emerald-800/20 border-emerald-300/40';
  if (id === 'obs' || id === 'var') return 'from-yellow-500/30 to-orange-700/20 border-amber-300/40';
  if (id === 'obsm' || id === 'varm') return 'from-sky-500/30 to-blue-700/20 border-sky-300/40';
  if (id === 'obsp' || id === 'varp') return 'from-violet-500/30 to-purple-800/20 border-violet-300/40';
  if (id === 'raw') return 'from-lime-500/30 to-emerald-700/20 border-lime-300/40';
  return 'from-slate-500/30 to-slate-700/20 border-slate-300/40';
};

export const DataStructureAtlas: React.FC<{ onBack: () => void }> = ({ onBack }) => {
  const [selectedId, setSelectedId] = useState<SlotId>('obs');
  const [query, setQuery] = useState('');

  const filtered = useMemo(() => {
    const q = query.trim().toLowerCase();
    if (!q) return SLOT_MAP;
    return SLOT_MAP.filter((item) =>
      [
        item.annPath,
        item.annLabel,
        item.seuPath,
        item.seuLabel,
        item.annMeaning,
        item.seuMeaning,
      ]
        .join(' ')
        .toLowerCase()
        .includes(q)
    );
  }, [query]);

  const selected = useMemo(
    () => filtered.find((item) => item.id === selectedId) ?? SLOT_MAP.find((item) => item.id === selectedId) ?? SLOT_MAP[0],
    [filtered, selectedId]
  );

  return (
    <div className="w-full h-screen overflow-y-auto bg-slate-950 text-white p-4 md:p-6">
      <div className="max-w-7xl mx-auto">
        <div className="flex flex-wrap items-center justify-between gap-3 mb-4">
          <button
            onClick={onBack}
            className="inline-flex items-center gap-2 rounded-lg border border-slate-700 bg-slate-900 px-3 py-2 text-sm hover:bg-slate-800"
          >
            <ArrowLeft className="w-4 h-4" />
            Back to Menu
          </button>
          <div className="inline-flex items-center gap-2 border border-slate-700 rounded-lg bg-slate-900 px-3 py-2 text-xs text-slate-300">
            <Sparkles className="w-4 h-4 text-cyan-300" />
            Visual object map for Seurat users learning AnnData
          </div>
        </div>

        <div className="rounded-xl border border-slate-800 bg-slate-900/70 p-4">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <h2 className="text-2xl font-bold flex items-center gap-2">
              <Database className="w-6 h-6 text-cyan-300" />
              AnnData ↔ Seurat v5 Structure Atlas
            </h2>
            <div className="text-xs text-slate-400">
              Click a block on either side to inspect mapping
            </div>
          </div>

          <div className="mt-3">
            <input
              type="text"
              value={query}
              onChange={(event) => setQuery(event.target.value)}
              placeholder='Search slots: "obs", "graphs", "layers", "umap"...'
              className="w-full rounded border border-slate-700 bg-slate-950 px-3 py-2 text-sm text-slate-100 placeholder:text-slate-500 focus:outline-none focus:border-cyan-500/70"
            />
          </div>

          <div className="mt-4 grid grid-cols-1 xl:grid-cols-2 gap-4">
            <div className="rounded-lg border border-cyan-500/30 bg-cyan-950/20 p-3">
              <div className="text-[11px] uppercase tracking-wider text-cyan-200 mb-2">AnnData object view</div>
              <div className="grid grid-cols-4 grid-rows-5 gap-2 min-h-[340px]">
                {SLOT_MAP.map((item) => {
                  const isSelected = selected.id === item.id;
                  const isFilteredOut = filtered.length > 0 && !filtered.some((f) => f.id === item.id);
                  return (
                    <button
                      key={`ann-${item.id}`}
                      onClick={() => setSelectedId(item.id)}
                      className={`rounded border bg-gradient-to-br px-2 py-1.5 text-left transition ${ANN_LAYOUT[item.id]} ${colorFor(item.id)} ${isSelected ? 'ring-2 ring-cyan-300/80' : ''} ${isFilteredOut ? 'opacity-30' : ''}`}
                      title={item.annPath}
                    >
                      <div className="text-[11px] font-mono text-slate-100">{item.annLabel}</div>
                    </button>
                  );
                })}
              </div>
            </div>

            <div className="rounded-lg border border-violet-500/30 bg-violet-950/20 p-3">
              <div className="text-[11px] uppercase tracking-wider text-violet-200 mb-2">Seurat object view</div>
              <div className="grid grid-cols-3 grid-rows-4 gap-2 min-h-[340px]">
                {SLOT_MAP.map((item) => {
                  const isSelected = selected.id === item.id;
                  const isFilteredOut = filtered.length > 0 && !filtered.some((f) => f.id === item.id);
                  return (
                    <button
                      key={`seu-${item.id}`}
                      onClick={() => setSelectedId(item.id)}
                      className={`rounded border bg-gradient-to-br px-2 py-1.5 text-left transition ${SEU_LAYOUT[item.id]} ${colorFor(item.id)} ${isSelected ? 'ring-2 ring-violet-300/80' : ''} ${isFilteredOut ? 'opacity-30' : ''}`}
                      title={item.seuPath}
                    >
                      <div className="text-[11px] font-mono text-slate-100 truncate">{item.seuLabel}</div>
                    </button>
                  );
                })}
              </div>
            </div>
          </div>

          <div className="mt-4 rounded-lg border border-slate-700 bg-slate-950/60 p-3">
            <div className="flex items-center gap-2 text-sm text-slate-200">
              <ArrowRightLeft className="w-4 h-4 text-cyan-300" />
              Selected mapping
            </div>
            <div className="mt-2 grid grid-cols-1 lg:grid-cols-2 gap-3">
              <div className="rounded border border-cyan-500/30 bg-cyan-950/20 p-3">
                <div className="text-[11px] uppercase tracking-wider text-cyan-200">AnnData slot</div>
                <div className="font-mono text-cyan-100 mt-1">{selected.annPath}</div>
                <p className="text-xs text-slate-200 mt-2">{selected.annMeaning}</p>
                <div className="mt-2 text-[11px] text-cyan-100 font-mono">{selected.shapeHint}</div>
                <pre className="mt-2 text-[11px] bg-slate-950/70 border border-slate-700 rounded px-2 py-1 overflow-x-auto text-slate-100">
{selected.pythonPeek}
                </pre>
              </div>
              <div className="rounded border border-violet-500/30 bg-violet-950/20 p-3">
                <div className="text-[11px] uppercase tracking-wider text-violet-200">Seurat location</div>
                <div className="font-mono text-violet-100 mt-1">{selected.seuPath}</div>
                <p className="text-xs text-slate-200 mt-2">{selected.seuMeaning}</p>
                <div className="mt-2 text-[11px] text-emerald-200">{selected.tip}</div>
                <pre className="mt-2 text-[11px] bg-slate-950/70 border border-slate-700 rounded px-2 py-1 overflow-x-auto text-slate-100">
{selected.rPeek}
                </pre>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};
