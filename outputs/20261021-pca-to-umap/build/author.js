const r = await artifactSession.run(async ({presentation})=>{
 const fs=await import('node:fs/promises');
 const root='/Users/ashi/github/informatics_club/outputs/20261021-pca-to-umap';
 const data=JSON.parse(await fs.readFile(root+'/build/numerical-example.json','utf8'));
 const slides=JSON.parse(await fs.readFile(root+'/build/revised-slides.json','utf8'));
 presentation.slides.keep([]);presentation.theme.defaultFont='Arial';
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
 // Slots 1–3 are replaced with the exact source slide parts after export.
 for(let i=1;i<=3;i++)base(i,'Preserved source slide');
 const sl4=base(4,'A cell begins as a vector of molecular counts');
 table(sl4,[['Cell',...data.genes,'Total UMIs'],...data.cells.map((c,i)=>[c,...data.counts[i].map(String),data.totals[i].toLocaleString('en-US')])],64,157,1152,295,23,[.7,1,1,1,1,1.3]);
 sl4.speakerNotes.text='The table now includes the original A, B and C values plus three additional cells and one additional gene. The six-cell, four-gene example lets us calculate a genuine four-to-three PCA reduction before a three-to-two UMAP embedding. All values are hypothetical, and total UMIs include omitted genes.\n\nجدول اکنون مقدارهای اولیهٔ A، B و C را همراه با سه سلول و یک ژن اضافه دارد. مثال شش‌سلولی و چهارژنی اجازه می‌دهد کاهش واقعی چهار به سه در PCA و سپس سه به دو در UMAP را محاسبه کنیم. همهٔ مقدارها فرضی‌اند و شمارش کل UMI شامل ژن‌های نشان‌داده‌نشده است.';
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
 await presentation.render({slide:[3,4,5,7,8,9,11,13,14,16,17,18,21],scale:1});
 await presentation.export({format:'pptx'});
 nodeRepl.write({slides:presentation.slides.items.length});
},{artifactType:'presentation',artifactRef:'artifact_HECR5XCkunqRaO7_Lu8T5MtdQ4Cf3S-R'});
nodeRepl.write({status:r.status,artifactRef:r.artifactRef,outputDirectory:r.outputDirectory,files:r.files});
