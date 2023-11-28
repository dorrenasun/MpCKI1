#################################### Environment ###############################
check_package<-function(pkg){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",update = T)
  if (!require(pkg,character.only = T,quietly = T)){
    # message(paste(pkg,"is not here. Install."))
    BiocManager::install(pkg,ask=F)
    library(pkg,character.only = T)
  }
}
pkg_list<-c("data.table","R.utils")
for (p in pkg_list) {check_package(p)}

#################################### Function ##################################
#Copy scripts
myCopy<-function(from,to,overwrite=recursive,recursive = FALSE,
                  copy.mode = TRUE, copy.date = FALSE){
  if (file.copy(from,to,overwrite = T)){
    message(paste(from,"copied to",to))
  }else{
    message(paste("Failed to copy",from))
  }
}

#Read all fasta files
read_fasta<-function(file){
  if (!file.exists(file)) stop("Fasta file not found!")
  Fasta<-fread(file,header = F,sep="\n",col.names = "Full")
  ID.positions<-grep(">",Fasta$Full,fixed = T)
  ID.lengths<-c(ID.positions[-1],nrow(Fasta)+1)-ID.positions
  IsIsoetes<-sum(grepl("Isoetes",Fasta$Full))
  if (IsIsoetes>0){
    Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
                ID=rep(sub(".*(Itaiw_v1_.*_[0-9]+_.*R.).*","\\1",Full[ID.positions]),ID.lengths))]
  }else{
    Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
                ID=rep(sub("[ \t].*","",Full[ID.positions]),ID.lengths))]
  }
  
  Fasta[,`:=`(ID=sub(">","",ID,fixed = T),
              Sequence=Full)]
  Fasta[ID.positions,Sequence:=""]
  Fasta.trans<-Fasta[,.(Sequence=paste(Sequence,collapse = "")),by=.(Num,ID)]
  Fasta.trans[,`:=`(ID=gsub("|","+",ID,fixed = T))]
  # Fasta.trans[,ID:=gsub("[ |\t]","_",ID)]
  #Check if the last character of sequence is *
  Fasta.trans[,`:=`(LastChar=substr(Sequence,nchar(Sequence),nchar(Sequence))
  )]
  Fasta.trans[LastChar=="*",Sequence:=substr(Sequence,1,nchar(Sequence)-1)]
  Fasta.trans[,Text:=paste(paste0(">",ID),Sequence,sep = "\n")]
  
  return(Fasta.trans[,.(ID,Sequence,Text)])
}

get_primary<-function(fasta){
  
  list.IDSrc<-list.files(Genomedir,"ID2Source.txt",full.names = T)
  IDSrc<-data.table()
  for (l in list.IDSrc){
    temp<-fread(l)
    IDSrc<-rbind(IDSrc,temp)
  }
  names(IDSrc)<-c("ID","Gene","var","IDSrc")
  fasta1<-copy(fasta)
  
  fasta1<-merge(fasta1,IDSrc,by="ID",all.x = T)
  #Make sure the variant is integer
  fasta1[,var:=as.integer(var)]
  fasta1[is.na(Gene),`:=`(Gene=ID,var=0)]
  #Convert variants to integers
  fasta1<-fasta1[order(Gene, var)]
  
  #Extract primary transcripts
  fasta.pri<-fasta1[!duplicated(Gene),]
  fasta.pri[,c("var","IDSrc"):=NULL]
  
  return(fasta.pri)
}

read_hmm<-function(hmmfile){
  hmm<-read.table(hmmfile,header = F,comment.char = "#")
  names(hmm)<-c("Target_name","Target_Acc","Query_name","Query_acc",
                "Evalue.Full","Score.Full","Bias.Full",
                "Evalue.Best1","Score.Best1","Bias.Best1",
                "Num_exp","Num_reg","Num_clu","Num_ov","Num_env","Num_dom","Num_rep","Num_inc",
                "Description"
  )
  return(as.data.table(hmm))
  
}
read_hmm_dom<-function(hmmdomfile){
  hmmdom<-read.table(hmmdomfile,header = F,comment.char = "#")
  names(hmmdom)<-c("Target_name","Target_Acc","Tlen","Query_name","Query_acc","Qlen",
                   "Evalue.Full","Score.Full","Bias.Full",
                   "No_Dom","Num_Dom","cEvalue.Dom","iEvalue.Dom","Score.Dom","Bias.Dom",
                   "Hmm.From","Hmm.To", "Ali.From","Ali.To","Env.From","Env.To",
                   "Accuracy","Description"
  )
  return(as.data.table(hmmdom))
  
}

#################################### Main ######################################
order.clade<-c("Chlorophytes",
               "Streptophyte algae",
               "Mosses",
               "Liverworts",
               "Hornworts",
               "Lycophytes",
               "Ferns",
               "Gymnosperms",
               "Angiosperms",
               "Basal Angiosperms",
               "Monocots",
               "Magnoliids",
               "Eudicots",
               "Unknown")
order.db<-c("Genome","OneKP","Dong","Manual")
order.tag<-c("Core","Neighbor","Outgroup")
Mac_home<-system("echo ~",intern = TRUE)
Homedir<-file.path(getwd(),"..")
# Homedir<-file.path(Mac_home,"Library/CloudStorage/OneDrive-Personal/Marchantia/Manuscripts/202304-MpCKI1/02-Github_submission/repos/00-Phylogenetics/CKI1")
codedir_ori<-file.path(Homedir,"00-Codes")
indir_ori<-file.path(Homedir,"01-Input")
outdir_ori<-file.path(Homedir,"02-Output")
koutdir_ori<-file.path(Homedir,"03-KeyOutput")

if (Round==1){
  para<-fread(file.path(Homedir,"02-Output","R1_timestamp.txt"),header = F)
}else if (Round==2){
  para<-fread(file.path(Homedir,"02-Output","R2_timestamp.txt"),header = F)
}else{
  stop("Round info not found!")
}

codedir<-file.path(outdir_ori,para[V1=="runID",V2],"00-Codes")
indir<-file.path(outdir_ori,para[V1=="runID",V2],"01-Input")
outdir<-file.path(outdir_ori,para[V1=="runID",V2],"02-Output")
koutdir<-file.path(outdir_ori,para[V1=="runID",V2],"03-KeyOutput")

if (!dir.exists(indir)) dir.create(indir)
if (!dir.exists(outdir)) dir.create(outdir)
if (!dir.exists(koutdir)) dir.create(koutdir)

OneKPdir<-dirname(para[V1=="Path_OneKP",V2])
Genomedir<-dirname(para[V1=="Path_Genome",V2])
Dongdir<-dirname(para[V1=="Path_Dong",V2])

# #Substitute the home directory in paths
# OneKPdir<-sub("/Users/.*/Library",file.path(Mac_home,"Library"),OneKPdir)
# Genomedir<-sub("/Users/.*/Library",file.path(Mac_home,"Library"),Genomedir)
# Dongdir<-sub("/Users/.*/Library",file.path(Mac_home,"Library"),Dongdir)

thisfile<-file.path(codedir_ori,"00-Functions.R")
myCopy(thisfile,codedir)
