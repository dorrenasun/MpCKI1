rm(list=ls())
message("Running 07-ProcessIPR.R !")
Round=2
source("00-Functions.R")

#################################### Main ######################################

#Input files
iprdir_ori<-file.path(koutdir_ori,"Interproscan")
iprdir<-file.path(koutdir,"Interproscan")
if (!dir.exists(iprdir_ori)){
  stop("IPR file not found in 03-KeyOutput! Exit!")
} else{
  if (!dir.exists(iprdir)) dir.create(iprdir)
  list.ipr<-list.files(path = iprdir_ori,full.names = T)
  for (f in list.ipr){
    myCopy(f,iprdir)
  }
  list.tsv<-grep("\\.tsv$",list.ipr,value = T)
  
}

file.entry<-file.path(indir_ori,"Interpro_entry.list")
if (!file.exists(file.entry)){
  stop("Entry file not found in 01-InOutput! Exit!")
}else{
  myCopy(file.entry,indir,overwrite = T)
}

file.shortname<-file.path(indir_ori,"short_names.dat")
if (!file.exists(file.shortname)){
  stop("Shortname file not found in 01-InOutput! Exit!")
}else{
  myCopy(file.shortname,indir,overwrite = T)
}

#Load raw results
data.all<-data.table()
for (f in list.tsv){
  data.tmp<-fread(f,header = F,sep="\t",fill = T)
  data.all<-rbind(data.all,data.tmp)
}

#Attach column names: https://interproscan-docs.readthedocs.io/en/latest/OutputFormats.html
names(data.all)<-c(
  "query",#Protein accession (e.g. P51587)
  "md5",# Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
  "len",# Sequence length (e.g. 3418)
  "db_info",# Analysis (e.g. Pfam / PRINTS / Gene3D)
  "db_acc",# Signature accession (e.g. PF09103 / G3DSA:2.40.50.140)
  "db_annot",# Signature description (e.g. BRCA2 repeat profile)
  "start",# Start location
  "stop",# Stop location
  "score",# Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
  "status",# Status - is the status of the match (T: true)
  "date",# Date - is the date of the run (dd-mm-yyyy)
  "ipr_acc",# InterPro annotations - accession (e.g. IPR002093)
  "ipr_annot",# InterPro annotations - description (e.g. BRCA2 repeat)
  "go",# (GO annotations (e.g. GO:0005515) - optional column; only displayed if –goterms option is switched on)
  "pathway"# (Pathways annotations (e.g. REACT_71) - optional column; only displayed if –pathways option is switched on)"Query","md5")
)

#Filtering: extract the TMMHelix annotations
tm<-unique(data.all[db_info=="TMHMM",])

#Filtering: keep only the domain IPR entries

entry<-fread(file.entry)
data.ipr<-data.all[ipr_acc %in% entry[ENTRY_TYPE=="Domain",ENTRY_AC],]

#Merge the region of same IPR entry
data.ipr[,ident:=1:.N]

data.ipr.long<-data.ipr[,.(query,len,ipr_acc,start,stop,on_aa=start:stop),by=ident]
data.ipr.aa<-unique(data.ipr.long[,.(query,len,ipr_acc,on_aa)])

findstarts<-function(numbers){
  num_sorted<-sort(numbers)
  num_prev<-c(-1,num_sorted[-length(num_sorted)])
  num_diff<-num_sorted-1-num_prev
  starts<-num_sorted[num_diff!=0]
  return(starts)
}
findends<-function(numbers){
  num_sorted<-sort(numbers)
  num_next<-c(num_sorted[-1],-1)
  num_diff<-num_sorted+1-num_next
  ends<-num_sorted[num_diff!=0]
  return(ends)
}

data.ipr.merge1<-data.ipr.aa[,.(start=findstarts(on_aa),stop=findends(on_aa)),by=.(query,len,ipr_acc)]
data.ipr.merge1[,dom_len:=stop-start+1]
#Check the overlap between different domains. If a short domain has >80% percent overlap with a long domain, then keep the long one
data.ipr.merge2<-data.table()
for (q in unique(data.ipr.merge1$query)){
  data.sub<-data.ipr.merge1[query==q,]
  setorder(data.sub, -dom_len,start)
  data.sub.new<-data.sub[1,]
  if (nrow(data.sub)>1){
    for (i in 2:nrow(data.sub)){
      N_overlap<-0
      for (j in 1:nrow(data.sub.new)){
        s1=data.sub[i,start]
        e1=data.sub[i,stop]
        s2=data.sub.new[j,start]
        e2=data.sub.new[j,stop]
        overlap<-length(intersect(s1:e1,s2:e2))/data.sub[i,dom_len]
        if (overlap>0.8){
          N_overlap<-N_overlap+1
        }
      }
      if (N_overlap==0) data.sub.new<-rbind(data.sub.new,data.sub[i,])
      
          }
    
  }
  data.ipr.merge2<-rbind(data.ipr.merge2,data.sub.new)
  
}

#Load short names for interpro
shortname<-fread(file.shortname,header = F)
names(shortname)<-c("ipr_acc","shortname")
data.ipr.merge3<-merge(data.ipr.merge2,shortname,by= "ipr_acc")
setorder(data.ipr.merge3,query,start)

#Final output
data.final<-rbind(tm[,.(query,len,db=db_info,acc=db_acc,name=db_acc,start,stop)],
                  data.ipr.merge3[,.(query,len,db="interpro",acc=ipr_acc,name=shortname,start,stop)]
                  )
setorder(data.final,query,start)
data.uniq<-data.final[!duplicated(acc),]

file.out<-file.path(outdir,"07-Subtree-IPR.tsv")
write.table(data.final,file.out,sep = "\t",row.names = F,col.names = T)
myCopy(file.out,koutdir_ori,overwrite = T)
