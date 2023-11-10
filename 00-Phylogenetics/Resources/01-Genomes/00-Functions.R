if (!require(data.table)) install.packages("data.table")
# if (!require(R.utils)) install.packages("R.utils")
library(data.table)
# library(R.utils)

#Read all fasta files
read_fasta<-function(file){
  if (!file.exists(file)) stop("Fasta file not found!")
  Fasta<-fread(file,header = F,sep="\t",col.names = "Full")
  ID.positions<-grep(">",Fasta$Full,fixed = T)
  ID.lengths<-c(ID.positions[-1],nrow(Fasta)+1)-ID.positions
  IsCoGe<-sum(grepl("||",Fasta$Full,fixed = T))
  if (IsCoGe>0){
    message(paste0(file," is a CoGe file! Parsing!"))
    Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
                       ID=rep(sub("^.*?\\|\\|.*?\\|\\|.*?\\|\\|.*?\\|\\|", "",Full[ID.positions]),ID.lengths))]
    Fasta[,ID:=sub("\\|.*\\|frame","_fr",ID)]
  }else{
    Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
                ID=rep(sub(" .*","",Full[ID.positions]),ID.lengths))]
  }
  
  Fasta[,`:=`(ID=sub(">","",ID,fixed = T),
              Sequence=Full)]
  Fasta[ID.positions,Sequence:=""]
  Fasta.trans<-Fasta[,.(Sequence=paste(Sequence,collapse = "")),by=.(Num,ID)]
  Fasta.trans[,`:=`(ID=gsub("|","+",ID,fixed = T)#,
                    #Sequence=gsub("*","",Sequence,fixed = T)
                    )]
  Fasta.trans[,ID:=sub("[ |\t].*","",ID)]
  Fasta.trans[,Text:=paste(paste0(">",ID),Sequence,sep = "\n")]
  #Remove duplicated sequences with the same name and sequence
  Fasta.uniq<-Fasta.trans[!duplicated(Text),]
  #Find IDs that are still duplicated
  list.dup<-Fasta.uniq[duplicated(ID),ID]
  if (length(list.dup)>0){
    message(paste("Some IDs are duplicated in",file,"! Postfixes added!"))
    setnames(Fasta.uniq,"ID","oldID")
    Fasta.uniq[oldID %in% list.dup,ID:=paste(oldID,1:.N,sep = "_"),by=oldID]
    Fasta.uniq[is.na(ID),ID:=oldID]
    Fasta.uniq[,Text:=paste(paste0(">",ID),Sequence,sep = "\n")]
    
  }
  
  return(Fasta.uniq[,.(Num,ID,Sequence,Text)])
}

get_primary<-function(fasta){
  
  list.IDSrc<-c(paste0(outdir,"01-Genomes-ID2Source.txt")
                )
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

