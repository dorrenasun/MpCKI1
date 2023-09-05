rm(list = ls())
#################### Environment ####################
source("00-Common.R")

if (!require(lasso2)) {
  packageurl <-"https://cran.r-project.org/src/contrib/Archive/lasso2/lasso2_1.2-22.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

pkg_list<-c("data.table","DESeq2","tximeta","readxl",#"DEGreport",
            "eulerr","ggplot2","scales","ggrepel","ggbreak")
for (p in pkg_list) {check_package(p)}
# if(!require(org.Mpolymorpha.eg.db)) source("00-GO_clustering.R")


#################Parameters#######################

outDir<-file.path(path.out)
if (dir.exists(outDir)) unlink(sub("\\/$","",outDir), recursive=TRUE)
if (!dir.exists(outDir)) dir.create(outDir,recursive = T)

SpecialGenes<-c(MpCKI1="Mp2g02750",
                MpCKX2="Mp5g10910",
                MpRRB="Mp4g20600",
                MpGCAM1="Mp6g04830",
                MpCPS="Mp2g07200",
                MpKS="Mp6g05950",
                MpKAOL1="Mp4g23680",
                MpKAOL2="Mp1g25410",
                MpTAA="Mp5g14320",
                MpTRIHELIX28="Mp2g16340",
                MpTRIHELIX39="Mp7g12160",
                MpKNOX1="Mp5g01600",
                MpAGO4b="Mp6g20400",
                
                MpBELL5="Mp5g11060",
                MpBZR3="Mp2g23000",
                MpBNB="Mp3g23300",
                # MpHYPNOS="Mp5g18910",
                MpEC="Mp5g18000"
                )


################# Functions #######################
myFisher<-function(A,B,v,total){
  message(paste0("Calculating Fishe's exact test of whether ", B, "is enriched of ",A))
  # A<-"Mpcki1_Veg.Up"
  # B<-"Mprrb_Veg.Up"
  total<-length(row.names(dds.All.kept))
  ori<-as.data.table(v$original.values,keep.rownames = T)
  names(ori)<-c("Name","Num")
  ori[,newName:=paste0(sort(unlist(strsplit(Name,"&"))),collapse = "&"),by="Name"]
  AnB<-paste0(sort(c(A,B)),collapse = "&")
  x<-ori[newName==AnB,Num]
  m<-ori[grepl(A,newName,fixed = T),sum(Num)]
  n<-total-m
  k<-ori[grepl(B,newName,fixed = T),sum(Num)]
  p<-phyper(x-1,m,n,k,lower.tail = FALSE)
  message(paste0("The p value is:",p))
  return(data.table(A=A,B=B,p=p))
}

myoutput<-function(DESeqResult,tag,SpecialGenes=NULL){
  file<-file.path(outDir,paste0("DEG-",tag,".tsv"))
  file.sig<-file.path(outDir,paste0("DEG-",tag,"-sig.tsv"))
  file.volcano<-file.path(outDir,paste0("Volcano-",tag,".pdf"))
  
  res <- suppressWarnings(as.data.table(DESeqResult,keep.rownames="Gene"))
  res[,Name:=Mpnames[match(Gene,MpTak_v6.1),MPGENES]]
  res[is.na(Name),Name:=""]
  if (!is.null(names(SpecialGenes))){
    res[Gene %in% SpecialGenes,Name:=names(SpecialGenes)[match(Gene,SpecialGenes)]]
    
  }
  cols<-c("Gene","Name",setdiff(names(res),c("Gene","Name")))
  
  res.out<-res[,..cols]
  write.table(res.out,
              file=file,
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t")
  res[,Sig:=((!is.na(padj)) & (padj<pthre) & (abs(log2FoldChange)>FCthre))]
  
  res.sig<-res[Sig==TRUE,..cols]
  # res.sig[,Sig:=NULL]
  setkey(res.sig,log2FoldChange)
  write.table(res.sig,
              file=file.sig,
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t")
  
  #Volcanoplot
  #labels
  ylim<-25
  xlim<-12
  labels<-data.frame(
    x=c(xlim*0.5,-xlim*0.5),
    y=c(ylim*0.8,ylim*0.8),
    text=c(paste0("Up:\n",
                  nrow(res[Sig==T&log2FoldChange>0,])," genes"
    ),
    
    paste0("Down:\n",
           nrow(res[Sig==T&log2FoldChange<0,])," genes"
    )
    )
  )
  #Label MpCPS as orange
  res[,Cat:=as.integer(padj<pthre)]
  res[Sig==T,Cat:=2]
  res[Sig==T&log2FoldChange<0&Gene%in%SpecialGenes,Cat:=3]
  res[Sig==T&log2FoldChange>0&Gene%in%SpecialGenes,Cat:=4]
  # res[,Cat:=as.factor(Cat)]
  
  #Fix all p>ylim to ylim
  res[-log10(padj)>ylim,padj:=10^(-ylim-1)]
  
  p.vol<-ggplot(res[!is.na(padj),],aes(x=log2FoldChange,y=-log10(padj)))+
    geom_point(data=res[Cat<=2,],aes(color=as.factor(Cat)),size=1*Lwd,alpha=0.6,stroke=0)+
    geom_point(data=res[Cat>2,],aes(color=as.factor(Cat)),size=1*Lwd,alpha=1,stroke=0)+
    geom_text(data=labels,aes(x,y,label=text),hjust=0.5,size=6*Lwd)+
    
    geom_text_repel(data=res[Cat>2,],
                    aes(x=log2FoldChange,
                        y=-log10(padj),
                        label=sub("\\[(.*)\\].*","\\1",Name),
                        color=as.factor(Cat)),
                    min.segment.length =0,
                    segment.size=0.2*Lwd,
                    max.overlaps=500,
                    # hjust=0.5,nudge_y = 1,
                    size=6*Lwd)+
    # geom_vline(xintercept = 0,size=0.5*Lwd,color="blue")+
    # geom_hline(yintercept = -log10(pthre),size=0.5*Lwd,color="blue")+
    geom_hline(yintercept = ylim+1,size=1*Lwd,color="grey50")+
    theme_pub+theme(aspect.ratio = 1.0,legend.position = "none")+
    scale_color_manual(values = c("#FFF5AC","grey70","black","#01549D","#A3071E"))+
    scale_x_continuous(limits = c(-xlim,xlim),oob = squish)+
    scale_y_continuous(expand = c(0,0),limits = c(-.5,ylim+2),breaks = seq(0,25,5),oob=squish)
  print(p.vol)
  ggsave(file.volcano,plot = p.vol,width = 5,height = 5,units = "cm")
  
  
  
  res.sigmat<-res[,.(Gene=Gene,
                     Up=(Sig==T&(log2FoldChange >FCthre)),
                     Down=(Sig==T&(log2FoldChange < -FCthre))
  )
  ]
  # setnames(res.sigmat,"Up",)
  names(res.sigmat)[2:3]<-paste(tag,names(res.sigmat)[2:3],sep=".")
  return(res.sigmat)
}

#Convert DESeq2 results to matrix
Res2List<-function(inFile,Tag="Data",rm=NULL){
  DEGdata<-fread(inFile,header = T,sep = "\t")
  Sig.up<-DEGdata[!is.na(padj) & padj<pthre & log2FoldChange>FCthre & (!Gene %in% rm),Gene]
  Sig.down<-DEGdata[!is.na(padj) & padj<pthre & log2FoldChange<(-FCthre) & (!Gene %in% rm),Gene]
  list.out<-list(Sig.up,Sig.down)
  names(list.out)<-paste0(Tag,c(".Up",".Down"))
  return(list.out)
}

################# Main #######################
message("Running DEG analysis for Aki 2022...")

{#Creat a sampleTable for Aki et al. 2022 samples
  path.Aki<-"./Salmon_fem_Aki"
  file.info.Aki<-file.path(path.Aki,"RunInfo.tsv")
  if (!file.exists(file.info.Aki)) source("00-Read_DDBJ.R")
  info.Aki<-fread(file.info.Aki)
  
  file.quant.Aki<-list.files(path.Aki,
                             pattern="\\.sf$",recursive = T,full.names = T)
  sampleTable.Aki<-data.table(files=file.quant.Aki,
                              acc=sub(".*(DRR[0-9]+).*","\\1",file.quant.Aki))
  sampleTable.Aki[,names:=info.Aki[match(acc,alias),sample_name]]
  sampleTable.Aki[,`:=`(Genotype=sub("_.*","",names),
                        Tissue="V")]
  
  #Factorize
  order.Aki<-c("Tak-1","CKX2-OX", "Tak-2",   "rrbKO" )
  sampleTable.Aki[,Genotype:=factor(Genotype,levels = order.Aki)]

  # rownames(sampleTable) <- sampleTable$names
}
{#Creat a sampleTable for Mpcki1 samples
  path.Bao<-"./Salmon_fem_Bao/"
  file.quant.Bao<-list.files(path=path.Bao,pattern="\\.sf$",recursive = T,full.names = T)
  
  #Creat a sampleTable
  sampleTable.Bao<-data.table(files=file.quant.Bao,
                              names=sub(".*((6|WT).*)\\/.*","\\1",file.quant.Bao)
  )
  sampleTable.Bao[,`:=`(Genotype=sub("-.*","",names),
                        Tissue=sub(".*-(R|V).","\\1",names)
  )]
  sampleTable.Bao[Genotype=="6",Genotype:="Mpcki1"]
  
  #Factorize
  order.Bao<-c("WT","Mpcki1")
  order.tis<-c("V","R")
  sampleTable.Bao[,`:=`(Genotype=factor(Genotype,levels = order.Bao),
                        Tissue=factor(Tissue,levels = order.tis)
  )]
  
  #sort by Tissue and Genotype, then combine them into Cond
  setorder(sampleTable.Bao,Tissue,Genotype)
  sampleTable.Bao[,Cond:=factor(paste0(Genotype,".",Tissue))]
}
{ #Merge two sample tables
  cols<-c("files","names","Genotype","Tissue")
  sampleTable.All<-rbind(sampleTable.Aki[,..cols],sampleTable.Bao[,..cols])
  sampleTable.All[,`:=`(Genotype=factor(Genotype,levels = c(order.Aki,order.Bao)),
                        Tissue=factor(Tissue,levels = order.tis)
                        )]
  sampleTable.All[,Cond:=factor(paste0(Genotype,".",Tissue))]
  
  #Determine colors for each genotype
  color.all<-c("black",pal[7],"grey50",pal[8],"grey50",pal[13])
  names(color.all)<-c(order.Aki,order.Bao)
  
}

{#Build DESeqDataSet & filter by total counts
  #Generate tx2gene list
  tx2gene.fem<-ids.v6[!grepl("MpVg",gene)]
  #Select samples
  sampleTable.select<-sampleTable.All[!Cond %in% c("Tak-1.V","CKX2-OX.V"),]
  se.All <- tximeta(sampleTable.select,type="salmon",skipMeta=T,txOut=F,tx2gene = tx2gene.fem)
  
  # Full interaction model
  dds.All <- DESeqDataSet(se.All, design = ~ Cond)
  dds.All.keep<-rowSums(counts(dds.All)) >= 10
  dds.All.kept <- dds.All[dds.All.keep,]
  
  #log-scale transformations
  vsd.All <- vst(dds.All.kept, blind=FALSE)
  # vsd_mat<-assay(vsd.All)
  # vsd_mat.dt<-as.data.table(vsd_mat,keep.rownames = "Gene")
  

  
  # data table for making UpSet plot
  UpSet.all<-data.table(Gene=unique(rownames(dds.All.kept)))

} 


{  #Visualization of sample distances
  pcaData <- plotPCA(vsd.All, intgroup=c("Genotype","Tissue","names"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  {
    p.pca<-ggplot(pcaData, aes(x=PC1,y=PC2, color=Genotype)) +
      geom_point(aes(shape=Tissue),size=4*Lwd,alpha=0.8) +
      geom_text(aes(label=names),size=2,show.legend = F,nudge_x = 1,hjust=0,vjust=0.5)+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      scale_color_manual(values = color.all)+
      scale_shape_manual(values =c(16,17)) +
      coord_fixed() + theme_pub + 
      theme(aspect.ratio = 1,
            legend.position = "right")
    print(p.pca)
    out.pdf<-file.path(outDir,"Sample_distance_Aki+Mpcki1_VST.pdf")
    ggsave(file=out.pdf, plot=p.pca, width=10, height=10,units = "cm")
  }
}

{ #Repeat the caculation of DEGs:
  
  # Ref:https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
  test.all <- DESeq(dds.All.kept)
  res.veg<-results(test.all, contrast = c("Cond","Mpcki1.V","WT.V"),alpha = pthre,independentFiltering = F)
  sigmat.veg<-myoutput(DESeqResult = res.veg,tag = "Mpcki1_Veg",SpecialGenes=SpecialGenes)
  
  res.rep<-results(test.all, contrast = c("Cond","Mpcki1.R","WT.R"),alpha = pthre,independentFiltering = F)
  sigmat.rep<-myoutput(DESeqResult = res.rep,tag = "Mpcki1_Rep",SpecialGenes=SpecialGenes)
  
  res.rrb<-results(test.all, contrast = c("Cond","rrbKO.V","Tak-2.V"),alpha = pthre,independentFiltering = F)
  sigmat.rrb<-myoutput(DESeqResult = res.rrb,tag = "Mprrb_Veg",SpecialGenes=SpecialGenes)
  
  
  UpSet.comp<-Reduce(function(x,y) merge(x = x, y = y, by = "Gene",all.x=T), 
                     list(UpSet.all,
                          sigmat.veg,
                          sigmat.rep,
                          sigmat.rrb
                     ))
  
  
}


{ #Count sums
  cols.all<-setdiff(names(UpSet.comp),c("Gene","Sum"))
  UpSet.comp[is.na(UpSet.comp)] <- F
  UpSet.comp.sum<-data.table(cond=cols.all,
                            sum=as.numeric(UpSet.comp[, lapply(.SD, sum),.SDcols=cols.all])
  )
  

}
#Venndiagram
# library(eulerr)
{ # euler plots
  list.comp<-list(
                 Mpcki1_Veg_vs_Rep=c("Mpcki1_Veg.Up","Mpcki1_Rep.Up","Mpcki1_Veg.Down","Mpcki1_Rep.Down"),
                 Mprrb_vs_Mpcki1=c("Mpcki1_Veg.Down","Mpcki1_Veg.Up","Mprrb_Veg.Down","Mprrb_Veg.Up")
                 
                 )
  
  for (l in names(list.comp)){ 
    cols.sel<-unlist(list.comp[l])
    data.venn<-as.data.frame(UpSet.comp[,..cols.sel])
    rownames(data.venn)<-UpSet.comp$Gene
    set.seed(123)

    v <- euler(data.venn)  
    labels<-UpSet.comp.sum[match(cols.sel, cond),paste0(cond,"\n(Total:",sum,")")]
    # error_plot(v)
    eulerr_options(edges=list(lwd=1),
                   labels=list(fontsize = 6,fontfamily="sans",font=3),
                   quantities=list(fontsize = 8),
                   padding=unit(2,"mm")
                   )
    p<-plot(v, 
            #fills = c("#ffd5d5","#d5e5ff","#ff8080","#87aade"),
            fills = F,
            labels=labels,
            
            quantities = TRUE)
    print(p)
    ggsave(file.path(outDir,paste0("Euler-",l,".pdf")),plot = p,width = 6,height = 6,units = "cm")
    
    #Try to find genes in each sector
    list.sect<-as.data.table(v$original.values,keep.rownames = T)
    names(list.sect)<-c("Sect","Num")
    setorder(list.sect,-Num)
    cts<-UpSet.comp[,c("Gene",cols.sel),with=F]
    
    for (s in list.sect[Num>0,Sect]){
      s.conds<-unlist(strsplit(s,"&"))
      s.otherconds<-setdiff(cols.sel,s.conds)
      cts[,ct1:=rowSums(.SD,na.rm=T),.SDcols=s.conds]
      if (length(s.otherconds)>0){
        cts[,ct2:=rowSums(.SD,na.rm=T),.SDcols=s.otherconds]
      }else{
        cts[,ct2:=0]
      }
      
      
      
      cts[ct1==length(s.conds)&ct2==0,Tag:=s]
      
      if (nrow(cts[ct1==length(s.conds)&ct2==0,])!=list.sect[Sect==s,Num]) warning("Something wrong!")
    }
    
    # list.euler<-
    list.euler<-merge(cts[!is.na(Tag),c("Tag","Gene")],
                      Mpnames[,.(MpTak_v6.1,Name,annots)],
                      by.x="Gene",by.y="MpTak_v6.1",all.x=T
    )
    list.euler[!is.na(Name)&!Name=="",GeneName:=paste0(Gene," [",Name,"]")]
    # list.euler[!is.na(GeneName),GeneName:=Gene]
    list.euler[is.na(list.euler)]<-""
    
    list.euler[,Tag:=factor(Tag,levels = list.sect$Sect)]
    setorder(list.euler,Tag,Gene)
    list.euler[,No.:=1:.N,by=.(Tag)]
    write.table(list.euler[,.(Tag,No.,Gene,Name,GeneName,annots)],
                file=file.path(outDir,paste0("Euler-",l,"_GeneList.tsv")),
                row.names = F,col.names = T,
                sep = "\t",quote = T)
    
    
  }
  #Calculate Fishe's exact test
  if (l=="Mprrb_vs_Mpcki1"){
    
    p.Up<-myFisher(A="Mprrb_Veg.Up",B="Mpcki1_Veg.Up",v = v,total = length(row.names(dds.All.kept)))
    p.Up.1<-myFisher(A="Mpcki1_Veg.Up",B="Mprrb_Veg.Up",v = v,total = length(row.names(dds.All.kept)))
    p.Down<-myFisher(A="Mprrb_Veg.Down",B="Mpcki1_Veg.Down",v = v,total = length(row.names(dds.All.kept)))
    p.Down.1<-myFisher(A="Mpcki1_Veg.Down",B="Mprrb_Veg.Down",v = v,total = length(row.names(dds.All.kept)))
    
    p.comp<-Reduce(function(x,y) rbind(x=x, y = y), 
                       list(p.Up,p.Up.1,p.Down,p.Down.1
                       ))
    
    write.table(p.comp,file=file.path(outDir,"FisherTest.tsv"),sep = "\t",quote = F,row.names = F,col.names = T)
    
  }
}


message("DESeq comparison finished!")



