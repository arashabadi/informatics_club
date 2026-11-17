import fs from 'node:fs/promises';
import {pathToFileURL} from 'node:url';
const {PresentationFile,FileBlob}=await import(pathToFileURL('/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs').href);
const presentation=await PresentationFile.importPptx(await FileBlob.load('/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap/build/source-preserved.pptx'));

 const root='/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap';
 const data=JSON.parse(await fs.readFile(root+'/build/numerical-example.json','utf8'));
 const slides=JSON.parse(await fs.readFile(root+'/build/revised-slides.json','utf8'));
 presentation.slides.keep([0,1,2,3]);
 if(presentation.slides.items.length!==4) throw Error('Protected slide count mismatch');
 const C={ink:'#102A43',body:'#243B53',teal:'#007F82',muted:'#627D98',pale:'#EFF5F7',white:'#FFFFFF'};
 function tx(sl,t,x,y,w,h,size=30,color=C.body,bold=false){let sh=sl.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{width:0,fill:'none'}});sh.text=t;sh.text.style={fontSize:size,typeface:'Arial',color,bold};return sh;}
 function base(n,title){const s=presentation.slides.add({width:1280,height:720,background:{fill:C.white}});tx(s,title,64,40,1152,98,41,C.ink,true);tx(s,String(n).padStart(2,'0'),1140,674,75,24,17,C.muted);return s;}
 function table(sl,values,x=64,y=158,w=1152,h=355,size=27,widths){const t=sl.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnTracks:(widths||values[0].map(()=>1)).map(value=>({mode:'fr',value})),values});
 t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({fill:C.pale,textStyle:{fontSize:size,typeface:'Arial',color:C.body},margins:{top:7,bottom:7,left:13,right:10},borders:{bottom:{color:'#FFFFFF',width:1}}});
 t.cells.block({row:0,column:0,rowCount:1,columnCount:values[0].length}).assign({fill:C.teal,textStyle:{fontSize:size-1,typeface:'Arial',bold:true,color:C.white}});return t;}
 function scatter(sl,coords,x,y,w,h,xTitle,yTitle){
 const colors=['#007F82','#007F82','#B94D53','#E28B2D','#795BA7','#386B9D'];
 const series=[0,2,3,4,5].map(i=>({name:i===0?'A/B':data.cells[i],xValues:[coords[i][0]],values:[coords[i][1]],fill:colors[i],marker:{symbol:'circle',size:13,fill:colors[i]},dataLabels:{showSeriesName:true,position:'top',textStyle:{fontSize:22,fill:colors[i]}}}));
 const chart=sl.charts.add('scatter',{position:{left:x,top:y,width:w,height:h},series,hasLegend:false,scatterOptions:{style:'marker'},chartFill:C.white,plotAreaFill:C.white,xAxis:{title:xTitle,textStyle:{fontSize:17,fill:C.muted},majorGridlines:null},yAxis:{title:yTitle,textStyle:{fontSize:17,fill:C.muted},majorGridlines:{fill:'#DEE7ED',width:1}}});return chart;
 }
 // Preserve the source slides directly. Only the existing slide-4 table changes.
 const sl4=presentation.slides.items[3];const tab4=sl4.tables.items[0];
 if(tab4.rowCount!==4||tab4.columnCount!==5) throw Error('Unexpected source table dimensions');
 tab4.insertRows({index:4,count:3});tab4.insertColumns({index:4,count:1});
 const newValues=[['Cell',...data.genes,'Total UMIs'],...data.cells.map((c,i)=>[c,...data.counts[i].map(String),data.totals[i].toLocaleString('en-US')])];
 for(let row=0;row<7;row++)for(let col=0;col<6;col++)tab4.cells.set(row,col,newValues[row][col]);
 tab4.columnTracks=[.7,1,1,1,1,1.3].map(value=>({mode:'fr',value}));
 for(let row=0;row<7;row++)tab4.rows[row].height=295/7;
 tab4.cells.block({row:0,column:0,rowCount:7,columnCount:6}).assign({fill:C.pale,textStyle:{fontSize:23,typeface:'Arial',color:C.body},margins:{top:5,bottom:5,left:13,right:10}});
 tab4.cells.block({row:0,column:0,rowCount:1,columnCount:6}).assign({fill:C.teal,textStyle:{fontSize:23,typeface:'Arial',bold:true,color:C.white}});
 for(let n=5;n<=24;n++){
  if(n===23){let s=base(23,'The paper distinguishes summaries from displays');s.images.add({path:'/Users/ashi/github/informatics_club/outputs/20261021-count-umap/presentations/assets/Fig4.png',alt:'Original Figure 4',fit:'contain',position:{left:64,top:145,width:1152,height:505}});continue;}
  const c=slides.find(s=>s.number===n);const sl=base(n,c.title);sl.speakerNotes.text=c.notes;
  if(c.kind==='text'){
   c.lines.forEach((line,i)=>tx(sl,line,75,170+i*87,1120,75,32));
   if(c.equation)tx(sl,c.equation,75,550,1120,70,29,C.teal,true);
  }else if(c.kind==='scatter'){
   scatter(sl,data.embedding,60,151,720,465,'UMAP1','UMAP2');
   table(sl,c.values,810,166,400,389,23,[.8,1.1,1.1]);
  }else{
   let h=355,size=27;
   if(c.values.length===8){h=395;size=26;}
   if(c.values.length===7){h=370;size=27;}
   if(c.values[0].length===7)size=24;
   if(n===19||n===22){h=430;size=25;}
   if(n===17){h=266;size=26;}
   if(n===14){h=245;size=28;}
   if(n===11){h=325;size=28;}
   if(n===10){table(sl,c.values,64,158,720,375,25,[.7,1,1,1,1.35]);scatter(sl,data.scores,805,170,408,365,'PC1','PC2');tx(sl,'PC3 also enters UMAP',824,550,385,50,23,C.teal,true);}
   else{
    let widths;
    if(n===19)widths=[.85,1.2,1.6];if(n===22)widths=[1.05,1.25,1.25];if(n===13)widths=[.8,1,1,1.3,1];
    table(sl,c.values,64,158,1152,h,size,widths);
   }
   if(c.equation)tx(sl,c.equation,74,158+h+22,1135,80,n===9?28:30,C.teal,true);
   if(n===17)tx(sl,'Strong edges favor higher q. Absent edges favor lower q.',75,480,1120,80,31,C.teal,true);
  }
  if(c.bottom)tx(sl,c.bottom,74,627,1135,42,21,C.muted);
 }

for(let i=0;i<presentation.slides.items.length;i++){
 const blob=await presentation.export({slide:presentation.slides.items[i],format:'png',scale:1});
 await fs.writeFile(root+'/qa/new-'+String(i+1).padStart(2,'0')+'.png',new Uint8Array(await blob.arrayBuffer()));
}
await (await PresentationFile.exportPptx(presentation)).save(root+'/build/authored.pptx');
console.log(JSON.stringify({slides:presentation.slides.items.length,retainedSourceSlides:4,onlyTableChangedOnSlide4:true}));
