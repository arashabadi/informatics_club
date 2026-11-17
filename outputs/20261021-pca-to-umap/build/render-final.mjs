import fs from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
const dir='/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';
const {PresentationFile,FileBlob}=await import(pathToFileURL(dir).href);
const p=await PresentationFile.importPptx(await FileBlob.load('/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap/build/candidate.pptx'));
for(const i of Array.from({length:p.slides.items.length},(_,i)=>i)){const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});await fs.writeFile(`/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap/qa/final-${String(i+1).padStart(2,'0')}.png`,new Uint8Array(await b.arrayBuffer()));}
console.log('Rendered final deck, '+p.slides.items.length+' slides.');
