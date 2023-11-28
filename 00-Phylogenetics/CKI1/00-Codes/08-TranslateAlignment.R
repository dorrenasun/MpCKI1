rm(list=ls())
message("Running 08-TranslateAlignment.R !")
Round=2
source("00-Functions.R")
pkg_list<-c("treeio")
for (p in pkg_list) {check_package(p)}


#Copy scripts
thisfile<-file.path(codedir_ori,"08-TranslateAlignment.R")
myCopy(thisfile,codedir,overwrite = T)

#Load information for species etc.
file.id<-file.path(outdir,"06-Subtree.id")
if (!file.exists(file.id)){
  stop("File for sequence ids not found in 02-Output!Exit!")
}

# file.ids<-list.files(outdir,"02-Candidates.*id$",full.names = T)
ids<-fread(file.id)
#ids_long<-ids[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=.(Num,Name,Clade,Species,Database)]
SDcols<-setdiff(names(ids),"Name")
ids_long<-ids[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=SDcols]

#read tree file if exist
file.tree<-list.files(path = koutdir_ori,pattern = "08.*ordered\\.newick",full.names = T)
if (sum(grepl("_np",file.tree))>0){
  file.tree<-file.tree[grepl("_np",file.tree)][1]
}
tree<-read.newick(file.tree)

order.tree<-data.table(No=1:length(tree$tip.label),Name=gsub("[\"\']","",tree$tip.label))
order.tree.sub<-merge(order.tree[,.(Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=No],
                      ids_long,
                      by="Name_sep",
                      all.x=T
)
setorder(order.tree.sub,"No")



#read fasta alignment file
list.files<-list.files(outdir,"*.fasta",full.names = T)
list.files<-list.files[!grepl("for_human",list.files,fixed = T)]

for (f in list.files){
  if (!file.exists(f)) stop("Target file not found!")
  if (grepl(".phy$",f)){
    phy<-fread(f,header=T,sep = " ")
    names(phy)<-c("Code","Seq_seg")
    phy[,No:=1:.N]
    fasta<-phy[,.(Sequence=paste0(Seq_seg,collapse = "")),by=.(Code)]
    fasta[,newName:=ids_long$Name[match(Code,ids_long$Num)]]
    outname<-file.path(OutDir,sub(".phy",".forhuman.fasta",basename(f)))
    
  }
  if (grepl(".fasta",f)){
    fasta_ori<-read_fasta(f)
    fasta<-merge(fasta_ori,ids,by="ID",all.x = T)
    setnames(fasta,"Name","newName")
    if (sum(is.na(fasta$newName))>0){
      message("The following ids failed to be retrieved!")
      message(paste(fasta$ID[is.na(newName)],"\n"))
      fasta[is.na(newName),newName:=ID]
    }
    outname<-file.path(outdir,sub(".fasta","_for_human.fasta",basename(f)))
  }
  fasta[,newText:=paste0(">",newName,"\n",Sequence)]
  
  order.out<-c(intersect(order.tree.sub$ID,fasta$ID),setdiff(fasta$ID,order.tree.sub$ID))
  fasta[,ID:=factor(ID,order.out)]
  setorder(fasta,ID)
  
  write(fasta$newText,outname)
  message(paste0(f, " is translated! Outputfile is ",outname))
}
