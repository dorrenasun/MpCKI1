#################################### Environment ###############################
rm(list = ls())
message("Running 02-MergeFasta.R !")

if (!require(stringr)) install.packages("stringr")
if (!"stringr" %in% .packages()) library(stringr)

Round=1
source("00-Functions.R")

#Copy scripts
thisfile<-file.path(codedir_ori,"02-MergeFasta.R")
myCopy(thisfile,codedir)

inFile.fasta.man<-file.path(indir_ori,"Manual_seqs.fasta")

myCopy(inFile.fasta.man,indir,overwrite = T)

#################################### Function ##################################

merge_ID<-function(names,by=";"){
  names_sep<-setdiff(unlist(strsplit(names,by,fixed = T)),c(NA,"NA"))
  if (length(names_sep)>0){
    namelist<-data.table(Name=names_sep,
                         Man=(names_sep %in% manuals.id),
                         Len=nchar(names_sep)
                         )
    setorder(namelist,-Man,Len,Name)
    
    newname<-paste(namelist$Name,collapse = by)
    return(newname)
    
  }else{
    return(as.character(NA))
  }
}
#################################### Main ######################################


#Read fasta files

# list.man<-c("../01-Input/Subtree_seqs-manual.fasta")
OneKP.hit<-file.path(outdir,"01-OneKP-blast_hits.fasta")
Genome.hit<-file.path(outdir,"01-AnnotatedGenomes-blast_hits.fasta")
Dong.hit<-file.path(outdir,"01-Dong-blast_hits.fasta")

list.hit<-c(Genome.hit,OneKP.hit,Dong.hit)
fasta<-data.table()
for (f in list.hit){
  message(paste0("Reading ",f))
  temp<-read_fasta(f)
  temp$Source<-f
  fasta<-rbind(fasta,temp)
}

#Retrieve Species information
inFile.ID2Src.Genome<-list.files(path = Genomedir,pattern = "ID2Source",full.names = T)
inFile.ID2Src.OneKP<-list.files(path = OneKPdir,pattern = "ID2Source",full.names = T)
inFile.ID2Src.Dong<-list.files(path = Dongdir,pattern = "ID2Source",full.names = T)

myCopy(inFile.ID2Src.Genome,indir)
myCopy(inFile.ID2Src.OneKP,indir)
myCopy(inFile.ID2Src.Dong,indir)
ID2Src.Genome<-fread(inFile.ID2Src.Genome,sep = "\t")
ID2Src.OneKP<-fread(inFile.ID2Src.OneKP,sep = "\t")
ID2Src.Dong<-fread(inFile.ID2Src.Dong,sep = "\t")
ID2Src.all<-rbind(ID2Src.Genome[,.(ID=pID,Clade,Species,Code,Database=order.db[1])],
                  ID2Src.OneKP[,.(ID=newID1,Clade,Species,Code,Database=order.db[2])],
                  ID2Src.Dong[,.(ID=pID,Clade="Liverworts",Species,Code,Database=order.db[3])]
)



#Retrieve source information for all sequences
list.missing<-setdiff(fasta$ID,ID2Src.all$ID)
if (length(list.missing)>0) {
  message("Some IDs could not be found in ID2Src!") 
  message(paste(list.missing,sep = "\n"))
} 

fasta.full<-merge(fasta[,.(ID,Sequence)],ID2Src.all,by="ID",all.x = T)
fasta.full[Clade %in% c("Basal Angiosperms","Magnoliids","Eudicots","Monocots"),Clade := "Angiosperms"]



#Substitute some Clade names



{#Introduce manual sequences
  manuals<-read_fasta(inFile.fasta.man)
  #Default species info:NA
  manuals[,`:=`(Species=as.character(NA),
                Clade=as.character(NA),
                Code=as.character(NA),
                Source=inFile.fasta.man,
                Database="Manual")]
  #Correct for ID if there is species info
  setnames(manuals,"ID","ID_ori")
  manuals[grepl("\\[|\\]",ID_ori),`:=`(ID=sub("\\[.*\\]","",ID_ori),
                                   Taxa=sub(".*\\[(.*)\\].*","\\1",ID_ori)
  )]
  manuals[!grepl("\\[|\\]",ID_ori),ID:=ID_ori]
  manuals[,Text:=paste0(">",ID,"\n",Sequence)]
  
  #Load species info from input
  manuals[!is.na(Taxa),`:=`(Species=sub(".*\\+","",Taxa),
                           Clade=sub("\\+.*","",Taxa)
  )]
  manuals[!is.na(Taxa),Species:=gsub("_"," ",Species)]
}

#If manual seqs already included
SDcols<-names(fasta.full)
fasta.all<-rbind(fasta.full,manuals[,..SDcols])

manuals.id<-unique(unlist(strsplit(manuals$ID,";",fixed = T)))

#Merge duplicated sequence
fasta.uni<-fasta.all[,.(Name=merge_ID(c(ID),by=";"),
                        Clade=merge_ID(c(Clade),by=";"),
                        Species=merge_ID(c(Species),by=";"),
                        Code=merge_ID(c(Code),by=";"),
                        Database=merge_ID(c(Database),by=";")
),by=Sequence]

fasta.uni[Clade %in% c("Basal Angiosperms","Magnoliids","Eudicots","Monocots"),Clade := "Angiosperms"]

#Check if there are seqs that has NA info
list.NA<-fasta.uni[is.na(Name) | is.na(Clade)|is.na(Species)|is.na(Database),]
if (nrow(list.NA)>0){
  message("Some sequences have unknown info!") 
  message(paste(list.NA$Name,sep = "\n"))
}

#Mark the fragmental sequences (length<300)
fasta.uni[nchar(Sequence)<300,`:=`(Frag=TRUE,
                                   Name=paste0(Name,".f"))]
fasta.uni[is.na(Frag),Frag:=FALSE]

#Assign numbers to replace the IDs
fasta.uni[,`:=`(Num=str_pad(1:nrow(fasta.uni), 5, pad = "0"))]
#Mark full-length and fragmental sequences
fasta.uni[Frag==TRUE,Num:=paste0(Num,"S")] #S for short
fasta.uni[Frag==FALSE,Num:=paste0(Num,"L")] #L for long

#Combine numbers with the sequence
fasta.uni[,Text:=paste0(">",Num,"\n",Sequence)]

#Write files
write(fasta.uni$Text,file.path(outdir,"02-Candidates-all.fasta"))
write(fasta.uni[!Frag==TRUE,Text],file.path(outdir,"02-Candidates-all-fl.fasta"))
write(fasta.uni[Frag==TRUE,Text],file.path(outdir,"02-Candidates-all-frag.fasta"))
write(fasta.uni[Database=="Genome",Text],file.path(outdir,"02-Candidates-Annotated.fasta"))

write.table(fasta.uni[,.(Num,Name,Clade,Species,Code,Database)],file.path(outdir,"02-Candidates-all.id"),col.names = T,row.names = F,quote = F,sep = "\t")

#List of species
Sources=unique(fasta.full[,.(Clade,Species,Code,Database)])
Sources[,`:=`(Clade=factor(Clade,levels = order.clade),
              Database=factor(Database,levels = order.db)
              )]
setorder(Sources,Clade,Database,Species)

write.table(Sources,file.path(outdir,"02-Sources.tsv"),col.names = T,row.names = F,quote = F,sep = "\t")

message("Fastas merged!")