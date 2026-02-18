import React, { Suspense, lazy, useEffect, useState } from 'react';
import { AppMode } from './types';
import { LayoutGrid, Binary, ArrowRight, Code2, Database } from 'lucide-react';

const Explorer3D = lazy(() =>
  import('./components/Explorer3D').then((module) => ({ default: module.Explorer3D }))
);
const FormulaLab = lazy(() =>
  import('./components/FormulaLab').then((module) => ({ default: module.FormulaLab }))
);
const PythonKernelInputs = lazy(() =>
  import('./components/PythonKernelInputs').then((module) => ({ default: module.PythonKernelInputs }))
);
const DataStructureAtlas = lazy(() =>
  import('./components/DataStructureAtlas').then((module) => ({ default: module.DataStructureAtlas }))
);

export default function App() {
  const [mode, setMode] = useState<AppMode>('MENU');

  useEffect(() => {
    const prevBodyOverflow = document.body.style.overflow;
    const prevHtmlOverflow = document.documentElement.style.overflow;
    const overflow = mode === 'EXPLORER' ? 'hidden' : 'auto';
    document.body.style.overflow = overflow;
    document.documentElement.style.overflow = overflow;
    return () => {
      document.body.style.overflow = prevBodyOverflow;
      document.documentElement.style.overflow = prevHtmlOverflow;
    };
  }, [mode]);

  if (mode !== 'MENU') {
      return (
        <Suspense fallback={
          <div className="w-full h-screen bg-slate-950 text-slate-200 flex items-center justify-center text-sm">
            Loading module…
          </div>
        }>
          {mode === 'EXPLORER' && <Explorer3D onBack={() => setMode('MENU')} />}
          {mode === 'FORMULAS' && <FormulaLab onBack={() => setMode('MENU')} />}
          {mode === 'PYTHON_INPUTS' && <PythonKernelInputs onBack={() => setMode('MENU')} />}
          {mode === 'STRUCTURE_ATLAS' && <DataStructureAtlas onBack={() => setMode('MENU')} />}
        </Suspense>
      );
  }

  return (
    <div className="w-full min-h-screen bg-slate-950 text-white flex flex-col items-center justify-center relative overflow-hidden">
        {/* Background Effects */}
        <div className="absolute top-0 left-0 w-full h-full bg-[radial-gradient(circle_at_50%_50%,_rgba(16,185,129,0.1),_transparent_70%)] pointer-events-none" />
        
        <div className="z-10 text-center max-w-6xl px-6">
            <h1 className="text-6xl font-black tracking-tighter mb-4 bg-clip-text text-transparent bg-gradient-to-r from-blue-400 via-emerald-400 to-purple-400">
                CellRank 2
            </h1>
            <p className="text-xl text-slate-400 mb-12 font-light">
                Interactive Learning Suite for Single-Cell Dynamics
            </p>

            <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-4 gap-6 w-full">
                {/* Card 1: 3D Explorer */}
                <button 
                    onClick={() => setMode('EXPLORER')}
                    className="group relative bg-slate-900 border border-slate-800 hover:border-blue-500 p-8 rounded-2xl transition-all duration-300 hover:shadow-[0_0_30px_rgba(59,130,246,0.2)] text-left flex flex-col h-64"
                >
                    <div className="bg-blue-500/10 w-12 h-12 rounded-lg flex items-center justify-center mb-6 group-hover:scale-110 transition-transform">
                        <LayoutGrid className="text-blue-400" size={24} />
                    </div>
                    <h2 className="text-2xl font-bold text-white mb-2">Explorer Studio</h2>
                    <p className="text-slate-400 text-sm mb-auto">
                        3D manifold + live formulas together. Tune kernel parameters and instantly see edges, matrices, and walks update.
                    </p>
                    <div className="flex items-center text-blue-400 font-bold text-sm mt-4 group-hover:translate-x-2 transition-transform">
                        Launch <ArrowRight size={16} className="ml-2" />
                    </div>
                </button>

                {/* Card 2: Formula Lab */}
                <button 
                    onClick={() => setMode('FORMULAS')}
                    className="group relative bg-slate-900 border border-slate-800 hover:border-emerald-500 p-8 rounded-2xl transition-all duration-300 hover:shadow-[0_0_30px_rgba(16,185,129,0.2)] text-left flex flex-col h-64"
                >
                    <div className="bg-emerald-500/10 w-12 h-12 rounded-lg flex items-center justify-center mb-6 group-hover:scale-110 transition-transform">
                        <Binary className="text-emerald-400" size={24} />
                    </div>
                    <h2 className="text-2xl font-bold text-white mb-2">Formula Lab</h2>
                    <p className="text-slate-400 text-sm mb-auto">
                        Dedicated deep-dive modules for kernels, GPCCA, and absorption math without the 3D scene.
                    </p>
                    <div className="flex items-center text-emerald-400 font-bold text-sm mt-4 group-hover:translate-x-2 transition-transform">
                        Start Calculating <ArrowRight size={16} className="ml-2" />
                    </div>
                </button>

                {/* Card 3: Python Inputs */}
                <button 
                    onClick={() => setMode('PYTHON_INPUTS')}
                    className="group relative bg-slate-900 border border-slate-800 hover:border-violet-500 p-8 rounded-2xl transition-all duration-300 hover:shadow-[0_0_30px_rgba(139,92,246,0.2)] text-left flex flex-col h-64"
                >
                    <div className="bg-violet-500/10 w-12 h-12 rounded-lg flex items-center justify-center mb-6 group-hover:scale-110 transition-transform">
                        <Code2 className="text-violet-400" size={24} />
                    </div>
                    <h2 className="text-2xl font-bold text-white mb-2">Python Input Arcade</h2>
                    <p className="text-slate-400 text-sm mb-auto">
                        Explore real AnnData structure snapshots (obs, layers, obsm, obsp) for every CellRank kernel.
                    </p>
                    <div className="flex items-center text-violet-400 font-bold text-sm mt-4 group-hover:translate-x-2 transition-transform">
                        Open Recipes <ArrowRight size={16} className="ml-2" />
                    </div>
                </button>

                {/* Card 4: Structure Atlas */}
                <button
                    onClick={() => setMode('STRUCTURE_ATLAS')}
                    className="group relative bg-slate-900 border border-slate-800 hover:border-cyan-500 p-8 rounded-2xl transition-all duration-300 hover:shadow-[0_0_30px_rgba(34,211,238,0.2)] text-left flex flex-col h-64"
                >
                    <div className="bg-cyan-500/10 w-12 h-12 rounded-lg flex items-center justify-center mb-6 group-hover:scale-110 transition-transform">
                        <Database className="text-cyan-300" size={24} />
                    </div>
                    <h2 className="text-2xl font-bold text-white mb-2">Structure Atlas</h2>
                    <p className="text-slate-400 text-sm mb-auto">
                        Separate visual map for AnnData vs Seurat v5 object structure with interactive slot mapping.
                    </p>
                    <div className="flex items-center text-cyan-300 font-bold text-sm mt-4 group-hover:translate-x-2 transition-transform">
                        Open Atlas <ArrowRight size={16} className="ml-2" />
                    </div>
                </button>

            </div>
        </div>

        <div className="absolute bottom-8 text-slate-600 text-xs">
            Built for Biologists • CellRank 2 Educational Demo
        </div>
    </div>
  );
}
