#################################### Environment ###############################
rm(list=ls())
message("Running 05-TranslateTree.R !")
args <- commandArgs(trailingOnly = TRUE)
# args<-2
Round=as.numeric(args[1])
source("00-Functions.R")

#Copy scripts
thisfile<-file.path(codedir_ori,"05-TranslateTree.R")
myCopy(thisfile,codedir,overwrite = T)

#################################### Function ##################################
myReplace<-function(entry){
  # entry<-data$Original
  temp1<-gsub("([\\(:\\);\\|/,]+)","@\\1@",entry)
  temp2<-data.table(Split1=unlist(strsplit(temp1,"@",fixed = T)))
  temp2[Split1 %in% dict$Num,Replace:=dict[match(Split1,Num),ID]]
  temp2[is.na(Replace),Replace:=Split1]
  #Check if any problem with splitting
  mimic<-paste0(temp2$Split1,collapse = "")
  if (!mimic==entry) message("Something wrong with splitting!")
  out<-paste0(temp2$Replace,collapse = "")
  return(out)
}



#################################### Main ######################################

list<-list.files(outdir,pattern="\\.(contree|treefile|support)$",full.names = T,recursive = F)

#Make a dictionary
if (Round==1){
  OTUList<-file.path(outdir,"02-Candidates-all.id")
}else if (Round==2){
  OTUList<-file.path(koutdir_ori,"02-Candidates-all.id")
  myCopy(OTUList,koutdir)
}


if (file.exists(OTUList)){
  dict<-fread(OTUList,header=T,colClasses = "character")
  setnames(dict,"Name","ID")
  testlen<-Inf
  # suffix<-paste0(rep("T",testlen-max(nchar(dict$Num))),collapse = "")
  dict[nchar(ID)>testlen,ID:=substr(ID,1,testlen)]
  dict[duplicated(ID),ID:=Num]
  # dict[,ID:=gsub(";","|",ID)]
  # dict[,ID:=gsub("[\\(\\)]","",ID)]
  # dict[,ID:=paste0(suffix,Num)]
  dict[,ID:=paste0("'",ID,"'")]
}

#Process each file
for (file in list){
  message(paste0("Processing ",file))
  
  data<-fread(file,header=F,quote="",sep = "\t",stringsAsFactors = F,colClasses = "character")
  names(data)<-"Original"
  data[,No:=1:.N]
  
  data[,Replace:=myReplace(Original),by=No]
  output<-paste0(file,".newick")
  
  for (i in data$No){
    if (i>1) output<-paste0(output,".",i+1)
    write.table(data[No==i,Replace],output,quote = F,row.names = F,col.names = F) 
  }
  
}

