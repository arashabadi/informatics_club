import React, { useState } from 'react';
import { AppMode } from './types';
import { Explorer3D } from './components/Explorer3D';
import { FormulaLab } from './components/FormulaLab';
import { PythonKernelInputs } from './components/PythonKernelInputs';
import { LayoutGrid, Binary, ArrowRight, Code2 } from 'lucide-react';

export default function App() {
  const [mode, setMode] = useState<AppMode>('MENU');

  if (mode === 'EXPLORER') {
      return <Explorer3D onBack={() => setMode('MENU')} />;
  }

  if (mode === 'FORMULAS') {
      return <FormulaLab onBack={() => setMode('MENU')} />;
  }

  if (mode === 'PYTHON_INPUTS') {
      return <PythonKernelInputs onBack={() => setMode('MENU')} />;
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

            <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-6 w-full">
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
                        Explore real AnnData structure snapshots (obs, layers, obsm, obsp) plus runnable Python for every CellRank kernel.
                    </p>
                    <div className="flex items-center text-violet-400 font-bold text-sm mt-4 group-hover:translate-x-2 transition-transform">
                        Open Recipes <ArrowRight size={16} className="ml-2" />
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
