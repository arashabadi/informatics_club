import fs from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
const SKILL='/Users/ashi/.codex/plugins/cache/openai-primary-runtime/presentations/26.1020.11127/skills/presentations';
const root='/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap';
const {finalizePresentation}=await import(pathToFileURL(SKILL+'/container_tools/artifact_tool_utils.mjs').href);
const tables=[4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,21,22,27,36];
const result=await finalizePresentation({
 workspaceDir:root,candidatePath:root+'/build/candidate.pptx',finalPath:root+'/output/PCA-to-UMAP_final.pptx',
 pythonExecutable:'/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3',
 integrityValidatorPath:SKILL+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:SKILL+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',...tables.flatMap(n=>['--require-native-table-slide',String(n)])],
 requiredNativeTableOwnerSlides:tables,requiredNativeChartOwnerSlides:[10,18],
 materializeLiteralChartWorkbooks:true,nativeChartTargetApplication:'powerpoint',
 verifyArtifactToolImport:true,receiptPath:root+'/build/validation-latest.json'
});
console.log(JSON.stringify({file:result.finalPath,package:result.packageIntegrity?.status,warnings:result.presentationLayout?.warnings,import:result.firstPartyImport?.passed}));
