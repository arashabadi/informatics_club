const fs = await import('node:fs/promises');
const root = '/Users/ashi/github/informatics_club/outputs/20261021-count-umap/presentations';
const content = JSON.parse(await fs.readFile(root+'/build/deck-content.json','utf8'));
const result = await artifactSession.run(async ({presentation,session}) => {
  presentation.theme.defaultFont = 'Arial';
  const C={navy:'#102A43',teal:'#007F82',ink:'#243B53',muted:'#627D98',pale:'#EFF5F7',white:'#FFFFFF',orange:'#D96920'};
  const text=(sl,t,x,y,w,h,size=30,color=C.ink,bold=false,font='Arial')=>{
    const sh=sl.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
    sh.text=t; sh.text.style={typeface:font,fontSize:size,color,bold};return sh;
  };
  const img=(sl,path,x,y,w,h)=>sl.images.add({path,alt:path.split('/').pop().replace('.png',''),fit:'contain',position:{left:x,top:y,width:w,height:h}});
  const footer=(sl,s)=>{
    let src=s.label?`Luecken & Theis 2019, ${s.label}, CC BY 4.0`:(s.kind==='pairimage'?'Official UMAP tutorial, synthetic data':'');
    if(src)text(sl,src,64,674,1060,24,16,C.muted);
    text(sl,String(s.number).padStart(2,'0')+(s.appendix?' / Appendix':''),1085,674,150,24,16,C.muted);
  };
  const list=(sl,lines,x,y,w,size=31,gap=75)=>lines.forEach((l,i)=>text(sl,l,x,y+i*gap,w,gap-10,size,C.ink));
  for(const s of content){
    const dark=['cover','divider'].includes(s.kind);
    const sl=presentation.slides.add({width:1280,height:720,background:{fill:dark?C.navy:C.white}});
    sl.speakerNotes.text=s.notes;
    if(s.kind==='cover'){
      text(sl,'scTalk 8',70,65,1100,40,25,'#76CECB',true);
      text(sl,s.title,70,175,1150,210,60,C.white,true);
      text(sl,s.subtitle,73,423,1090,70,36,'#76CECB');
      text(sl,s.footer,73,579,1000,85,23,'#D9E5EC');
      continue;
    }
    if(s.kind==='divider'){
      text(sl,'Backup material',70,140,1140,70,27,'#76CECB');
      text(sl,s.title,70,247,1140,150,55,C.white,true);
      text(sl,s.subtitle,73,452,1100,120,31,'#D9E5EC');
      continue;
    }
    const shortTitles={
      2:'Preprocessing prepares cells for biological inference',
      6:'HVG selection and scaling define the gene space',
      8:'Many genes share fewer patterns of variation',
      14:'Neighbor search defines local cell relationships',
      15:'UMAP fits a layout to a weighted graph',
      19:'Clustering and UMAP use related representations',
      23:'Discussion: which decision changes the conclusion?',
      27:'Trajectory models add assumptions beyond UMAP',
      30:'Random data can produce apparent marker evidence',
      31:'Each processing stage has different applications',
      32:'The classical Seurat sequence',
      33:'The classical Scanpy sequence'
    };
    text(sl,shortTitles[s.number]||s.title,64,44,1155,91,42,C.navy,true);
    footer(sl,s);
    if(s.kind==='figure'){
      if(s.wide){img(sl,s.image,64,143,1152,504);}
      else{
        img(sl,s.image,64,145,770,500);
        list(sl,s.points||[],870,193,343,29,127);
      }
    }else if(s.kind==='wideimage'){
      img(sl,s.image,64,140,1152,446);
      text(sl,s.caption,73,607,1130,50,26,C.teal);
    }else if(s.kind==='pairimage'){
      img(sl,s.images[0],64,148,556,365);img(sl,s.images[1],660,148,556,365);
      text(sl,s.caption,64,528,1150,40,23,C.muted);
      list(sl,s.lines,64,578,1140,26,38);
    }else if(s.kind==='table'){
      const n=s.values.length;
      const height=n>=6?425:n===5?340:n===4?295:235;
      const table=sl.tables.add({rows:n,columns:s.values[0].length,left:64,top:157,width:1152,height,
        columnTracks:(s.widths||s.values[0].map(()=>1)).map(value=>({mode:'fr',value})),values:s.values});
      const all=table.cells.block({row:0,column:0,rowCount:n,columnCount:s.values[0].length});
      all.assign({fill:C.pale,textStyle:{fontSize:s.number===31?23:26,typeface:'Arial',color:C.ink},margins:{left:16,right:14,top:12,bottom:12},borders:{bottom:{width:1,color:'#FFFFFF'}}});
      table.cells.block({row:0,column:0,rowCount:1,columnCount:s.values[0].length}).assign({fill:C.teal,textStyle:{fontSize:25,typeface:'Arial',bold:true,color:C.white}});
      let y=157+height+27;
      if(s.formula){text(sl,s.formula,70,y,1130,52,31,C.teal);y+=63;}
      if(s.lines)list(sl,s.lines,70,y,1130,26,43);
      if(s.note)text(sl,s.note,70,629,1130,37,18,C.muted);
    }else if(s.kind==='code'){
      text(sl,s.code,70,151,1140,440,s.number===33?27:29,C.ink,false,'Menlo');
      text(sl,s.note,70,606,1135,59,22,C.muted);
    }else if(s.kind==='discussion'){
      list(sl,s.lines,70,177,1110,34,112);
      text(sl,s.note,70,565,1110,85,29,C.teal,true);
    }else if(s.kind==='sources'){
      list(sl,s.lines,70,146,1130,29,57);
      text(sl,s.small,70,404,1130,225,23,C.teal);
    }else{
      let y=161;
      if(s.formula){text(sl,s.formula,70,157,1130,76,44,C.teal,true);y=265;}
      const gap=s.lines.length===5?68:75;
      list(sl,s.lines,70,y,1130,31,gap);
      if(s.big)text(sl,s.big,70,s.formula?534:(s.lines.length===5?541:535),1130,86,s.number===8?32:31,C.teal,true);
      if(s.note)text(sl,s.note,70,634,1130,32,18,C.muted);
    }
  }
  nodeRepl.write({slides:presentation.slides.items.length,notes:content.length});
  const preview=await presentation.render({slide:[0,2,8,13,15,17,22,25,30,32],scale:1});
  nodeRepl.write(preview.map(x=>({name:x.name,size:x.size})));
}, {artifactType:'presentation',artifactRef:'artifact_HECR5XCkunqRaO7_Lu8T5MtdQ4Cf3S-R'});
nodeRepl.write(result);
