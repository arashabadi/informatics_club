import fs from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
const dir='/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';
const {PresentationFile,FileBlob}=await import(pathToFileURL(dir).href);
const p=await PresentationFile.importPptx(await FileBlob.load('/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap/build/source-preserved.pptx'));
for(const i of [0,1,2,3]){const b=await p.export({slide:p.slides.items[i],format:'png',scale:1});await fs.writeFile(`/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap/qa/source-${i+1}.png`,new Uint8Array(await b.arrayBuffer()));}
console.log('Rendered source slides 1–4.');
