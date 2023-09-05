#################################### Environment ###############################
rm(list=ls())
message("Running 06-ExtractFromSubtree.R !")
Round=2
source("00-Functions.R")
#Copy scripts
thisfile<-file.path(codedir_ori,"06-ExtractFromSubtree.R")
myCopy(thisfile,codedir,overwrite = T)

#################################### Function ##################################
read_tree<-function(nexfile){
  data<-fread(nexfile,header=F,quote="",sep = "\t",fill = T)
  Line_begin<-grep("begin trees;",data$V1)
  Line_end<-grep("end;",data$V1)
  trees<-data$V1[(Line_begin+1):(Line_end-1)]
  all_seqs<-data.table()
  
  for (entry in trees){
    #Remove all annotations from FigTree
    clean1<-gsub("\\[.*?\\]","",entry)
    clean2<-sub(".*=","",clean1)
    #Remove the scores and quotes
    clean3<-gsub(":[0-9\\.E\\-]+","",clean2)
    #Remove ending semicolon
    if (substr(clean3, nchar(clean3), nchar(clean3))==";"){
      clean3<-substr(clean3, 1, nchar(clean3)-1)
    }
    seqs<-data.table(Name_comb=unlist(strsplit(clean3,",")))
    #Remove leading and trailing spaces
    seqs[,Name_comb:=trimws(Name_comb,"both")]
    #Remove leading and trailing brackets
    seqs[,Name_comb:=trimws(Name_comb,"both",whitespace = "[\\(\\)]")]
    #Remove leading and trailing quotes
    seqs[,Name_comb:=trimws(Name_comb,"both",whitespace = "'")]
    
    seqs[,No:=1:.N]
    #Split the IDs by semicolons
    seqs_sep<-seqs[,.(No,Name_sep=unlist(strsplit(Name_comb,"[;\\|]"))),by=Name_comb]
    #Remove quotation marks
    seqs_sep[,Name_sep:=gsub("\\'","",Name_sep)]
    setkeyv(seqs,"No")
    all_seqs<-rbind(all_seqs,seqs_sep)
  }
  return(all_seqs)
}
#################################### Main ######################################
#Load files for subtrees
file.core<-list.files(path = koutdir_ori,pattern = "06-Subtree-Core.nex",full.names = T)
if (!file.exists(file.core)) stop("Subtree-Core file not found! Exit!")
myCopy(file.core,koutdir,overwrite = T)
core<-read_tree(file.core)
core[,ID:=paste0("C",No)]
core[,Tag:=order.tag[1]]

file.nei<-list.files(path = koutdir_ori,pattern = "06-Subtree-Neighbors.nex",full.names = T)
if (!file.exists(file.nei)) stop("Subtree-Neighbors file not found! Exit!")
myCopy(file.nei,koutdir,overwrite = T)
nei<-read_tree(file.nei)
nei[,ID:=paste0("N",No)]
nei[,Tag:=order.tag[2]]

list.iscore<-nei[Name_sep %in% core$Name_sep,Name_comb]
nei.notcore<-nei[!Name_comb %in% list.iscore,]

file.out<-list.files(path = koutdir_ori,pattern = "06-Subtree-Outgroup.nex",full.names = T)
if (!file.exists(file.out)) stop("Subtree-Uutgroup file not found! Exit!")
myCopy(file.out,koutdir,overwrite = T)
outgroup<-read_tree(file.out)
outgroup[,Tag:=order.tag[3]]
outgroup[,ID:=paste0("O",No)]

interested<-do.call("rbind",list(core,nei.notcore,outgroup))

#Manually remove sequences:
file.exclude<-list.files(path=indir_ori,pattern="Exclude.txt",full.names = T)
if (!file.exists(file.exclude)) {
  message("No file found for manual exclusion")
  exclude<-c()
  
}else{
  myCopy(file.exclude,indir,overwrite = T)
  exclude<-suppressWarnings(fread(file.exclude,header = F,sep = "\t",blank.lines.skip = T))
  if (nrow(exclude)>0){
    exclude<-setdiff(exclude[!grepl("#",V1),V1],"")
  }else {
      exclude<-c()
    } 
  #Check if these names are meaningful
  list.missing<-setdiff(exclude,interested$Name_comb)
  if (length(list.missing)>0){
    message("\nSome sequences to exlude not found:")
    message(paste(list.missing,collapse = "\n"))
  }
}

interested<-interested[!Name_comb %in% exclude,]

#Load file for general info
file.id<-list.files(path = koutdir_ori,pattern = "02-Candidates-all.id",full.names = T)
if (!file.exists(file.id)) stop("Fasta id file not found! Exit!")
myCopy(file.id,koutdir,overwrite = T)
id<-fread(file.id)

id[Database=="Genome;Manual",Database:="Genome"]
id.sub1<-id[,.(Database,Clade,Species_sep=unlist(strsplit(Species,"[;\\|]"))),by=.(Num,Species)]
id.sub2<-id[,.(Name_sep=unlist(strsplit(Name,"[;\\|]")),Num),by=Name]

#Find out if IDs should be included by species
file.include<-list.files(path = koutdir_ori,pattern = "02-Sources.tsv",full.names = T)
if (!file.exists(file.include)) stop("Source file not found! Exit!")
myCopy(file.include,koutdir,overwrite = T)
include<-fread(file.include)
include.id<-merge(id.sub1,include,by.x=c("Clade","Species_sep","Database"),by.y=c("Clade","Species","Database"),all.x = T,allow.cartesian = T)

#Find out IDs for names of interest
interested.id<-merge(interested,id.sub2[,.(Name_sep,Num)],by="Name_sep")
#Check
list.missing<-setdiff(interested$ID,interested.id$ID)
if (length(list.missing)>0) warning("Not all IDs were found!")

#Decide whether to keep each seq
keep.core<-include.id[!is.na(Include_in_Core) & Include_in_Core=="T",Num]
interested.id[Tag==order.tag[1] & Num %in% keep.core, Keep:=T]
keep.nei<-include.id[!is.na(Include_in_Neighbors) & Include_in_Neighbors=="T",Num]
interested.id[Tag==order.tag[2] & Num %in% keep.nei, Keep:=T]
keep.out<-include.id[!is.na(Include_in_Outgroup) & Include_in_Outgroup=="T",Num]
interested.id[Tag==order.tag[3] & Num %in% keep.out, Keep:=T]

kept<-unique(interested.id[Keep==T,.(Num,Name_comb,Tag)])





#Read the source fasta
file.fasta<-list.files(path = koutdir_ori,pattern = "02-Candidates-all.fasta",full.names = T)
if (!file.exists(file.fasta)) stop("Fasta file not found! Exit!")
myCopy(file.fasta,koutdir,overwrite = T)
fasta<-read_fasta(file.fasta)

fasta.human<-merge(fasta,id,by.x="ID",by.y="Num",all.x=T)

#Check if all IDs could be found
list.missing<-setdiff(kept$Num,fasta.human$ID)
if (length(list.missing)>0) warning("Not all IDs were found in fasta!")

#Select sequence & merge info
fasta.kept<-fasta.human[ID %in% kept$Num,]
fasta.final<-merge(fasta.kept,kept[,.(ID=Num,Tag)],by="ID")
fasta.final[,HumanText:=paste0(">",Name,"\n",Sequence)]
fasta.final[,`:=`(Database=factor(Database,order.db),
                  Clade=factor(Clade,order.clade),
                  Tag=factor(Tag,order.tag)
                  )]
setorder(fasta.final,ID)
write(fasta.final$Text,file.path(outdir,"06-Subtree.fasta"))

setorder(fasta.final,-Tag,Clade,Database,Species,Name)
write(fasta.final$HumanText,file.path(outdir,"06-Subtree_for_human.fasta"))
write.table(fasta.final[,.(ID,Name,Species,Clade,Database,Tag)],
            file.path(outdir,"06-Subtree.id"),
            col.names = T,row.names = F,sep = "\t",quote = F
            )

#Split files for InterproScan
MaxSeg<-100
nseq<-nrow(fasta.final)
fasta.final[,Count:=1:nseq]
for (i in 1:ceiling(nseq/MaxSeg)){
  start<-(i-1)*MaxSeg+1
  end<-min(nseq,MaxSeg*i)
  fasta.seg<-fasta.final[Count>=start & Count<=end,]
  outfile.seg<-file.path(outdir,paste0("06-Subtree-seg",i,".fasta"))
  write(fasta.seg$Text,outfile.seg)
}

message("Sequences extracted from subtree!")
