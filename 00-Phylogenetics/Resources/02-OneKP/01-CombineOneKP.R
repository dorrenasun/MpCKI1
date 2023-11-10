rm(list = ls())
if (!require(data.table)) install.packages("data.table")
if (!require(R.utils)) install.packages("R.utils")
if (!require(stringr)) install.packages("stringr")
if(!require(readxl)) install.packages("readxl")

library(data.table)
library(R.utils)
library(stringr)
library(readxl)

#Function to read all fasta files
read_fasta<-function(file){
  if (!file.exists(file)) stop("Fasta file not found!")
  Fasta<-fread(file,header = F,sep="\t",col.names = "Full")
  ID.positions<-grep(">",Fasta$Full,fixed = T)
  ID.lengths<-c(ID.positions[-1],nrow(Fasta)+1)-ID.positions
  
  Fasta[,`:=`(Num=rep(1:length(ID.positions),ID.lengths),
              ID=rep(sub(" .*","",Full[ID.positions]),ID.lengths))]
  Fasta[,`:=`(ID=sub(">","",ID,fixed = T),
              Sequence=Full)]
  Fasta[ID.positions,Sequence:=""]
  Fasta.trans<-Fasta[,.(Sequence=paste(Sequence,collapse = "")),by=.(Num,ID)]
  Fasta.trans[,`:=`(ID=gsub("|","+",ID,fixed = T))]
  Fasta.trans[,ID:=sub("[ |\t].*","",ID)]
  Fasta.trans[,ID:=sub("scaffold-","",ID,fixed = T)]
  Fasta.trans[,ID:=gsub("\\(","{",ID)]
  Fasta.trans[,ID:=gsub("\\)","}",ID)]
  Fasta.trans[grepl("PIUF",ID,fixed=T),ID:=sub("Pellia_sp._{cf_epiphylla_{L.}_Corda}","Pellia_cf._epiphylla",ID,fixed = T)]
  
  #Only keep the first two or three words in taxon name
  Fasta.trans[grepl("cf",ID,fixed = T),newID:=sub("^([A-Z]+\\-[0-9]+\\-[[:alpha:]]+_cf\\._[[:alpha:]]+).*","\\1",ID)]
  Fasta.trans[!grepl("cf",ID,fixed = T),newID:=sub("^([A-Z]+\\-[0-9]+\\-[[:alpha:]]+_[[:alpha:]]+).*","\\1",ID)]
  
  #Fasta.trans[,Text:=paste(paste0(">",newID),Sequence,sep = "\n")]
  
  return(Fasta.trans)
}

#List files downloaded
inDir<-"./Download"
list.all<-data.table(File=list.files(inDir,pattern = "translated-protein",full.names = T,recursive = T))
list.all[,Code:=sub(".*\\/([A-Z]+)-.*","\\1",File)]

#Retrieve species info
Info<-fread(file.path(inDir,"Sample-List-with-Taxonomy.tsv.csv"))
file.prefix<-"Prefixes.xlsx"
if (!file.exists(file.prefix)) stop(paste(file.prefix," not found! Stop!"))
prefixes<-as.data.table(read_excel(file.prefix,sheet = 1))
Info.full<-merge(Info,prefixes[,.(Order,Clade_new,Prefix)],by="Order")

#Select files to merge
Info.sele<-Info.full[!grepl("#",Clade_new,fixed = T),]
Info.sele[,File:=list.all[match(Info.sele$`1kP_Sample`,Code),File]]
setnames(Info.sele,"1kP_Sample","Code")
if (nrow(Info.sele[is.na(File),])>0) warning("Some files not found!")

#Merge selected files
list.file<-Info.sele$File
fasta.all<-data.table()
for (f in list.file){
  print(paste0("Processing: ",f))
  temp<-read_fasta(f)
  temp$Source<-f
  #Remove "scaffold" from gene names
  fasta.all<-rbind(fasta.all,temp)
}

#Check if there are redundant IDs
message("Checking redundancy in IDs..")
if (sum(duplicated(fasta.all$ID))>0) warning("Duplicated IDs! Wrong parsing!")

#Add clade prefixes to IDs
fasta.all[,Code:=sub("-.*","",ID)]
fasta.all<-merge(fasta.all,Info.sele[,.(Code,Prefix,Clade_new,Species)],all.x = T)
if (sum(is.na(fasta.all$Prefix)>0)) message("Sequences of unknown source!")
fasta.all[,newID1:=paste0(Prefix,"_",newID)]

#Confirm ID format with sampling
fasta.sample<-fasta.all[!duplicated(Source),]

#Truncate IDs longer than 50
strlim<-50
fasta.all[nchar(newID1)>strlim,newID1:=str_trunc(newID1,strlim,ellipsis = "..")]
#Check if there are still IDs longer than 50
if (max(nchar(fasta.all$newID1))>strlim) warning("IDs longer than 50!")
#Put the new ID together with the sequence
fasta.all[,Text:=paste(paste0(">",newID1),Sequence,sep = "\n")]

#Extract sample info from the final fasta
Taxa<-fasta.all[!duplicated(Code),.(Code,Clade=Clade_new,Species,Source)]

#Write files
{
  outDir<-"Extracted/"
  outbase<-unlist(strsplit(inDir,"[\\.\\/]+"))
  outbase<-"OneKP-Nonseed"
  
  list.old<-list.files(outDir,full.names = T,recursive = T)
  if (sum(grep(outbase,list.old))>0) file.remove(list.old[grep(outbase,list.old)])
  
  write(fasta.all$Text,paste0(outDir,outbase,"-merged.fasta"))
  write.table(fasta.all[,.(newID1,Source,Code,Clade=Clade_new,Species)],paste0(outDir,outbase,"-ID2Source.txt"),sep = "\t",quote = F,row.names = F,col.names = T)
  write.table(Taxa,paste0(outDir,outbase,"-Samples.tsv"),sep = "\t",quote = F,row.names = F,col.names = T)
  
}
message("OneKP Transcriptomes merged!")