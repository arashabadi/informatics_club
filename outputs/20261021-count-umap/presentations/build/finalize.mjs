import fs from 'node:fs/promises';
import { pathToFileURL } from 'node:url';
const SKILL='/Users/ashi/.codex/plugins/cache/openai-primary-runtime/presentations/26.1015.12040/skills/presentations';
const root='/Users/ashi/github/informatics_club/outputs/20261021-count-umap/presentations';
const { finalizePresentation }=await import(pathToFileURL(SKILL+'/container_tools/artifact_tool_utils.mjs').href);
await fs.mkdir(root+'/output',{recursive:true});
const result=await finalizePresentation({
  workspaceDir:root,
  candidatePath:root+'/build/candidate.pptx',
  finalPath:root+'/output/From-Count-Matrix-to-UMAP_2026-10-21_v3.pptx',
  pythonExecutable:'/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3',
  integrityValidatorPath:SKILL+'/container_tools/inspect_presentation_package_integrity.py',
  layoutValidatorPath:SKILL+'/container_tools/inspect_presentation_layout_geometry.py',
  layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...[3,5,6,14,18,19,21,31].flatMap(n=>['--require-native-table-slide',String(n)])],
  requiredNativeTableOwnerSlides:[3,5,6,14,18,19,21,31],
  fontPolicy:{basis:'design',families:['Arial','Menlo'],scriptFonts:{cs:'Arial Unicode MS'}},
  verifyArtifactToolImport:true,
  receiptPath:root+'/build/final-validation-v3.json'
});
console.log(JSON.stringify(result,null,2));
