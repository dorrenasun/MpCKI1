#################################### Environment ###############################

rm(list=ls())
message("Running 03-FilterbyHmm.R !")
Round=1
source("00-Functions.R")

#Copy scripts
thisfile<-file.path(codedir_ori,"03-FilterbyHmm.R")
myCopy(thisfile,codedir,overwrite = T)

#################################### Main ######################################

#Load inputs
fastafile<-file.path(outdir,"02-Candidates-all-fl.fasta")
fasta<-read_fasta(fastafile)
hmmdomfile<-file.path(outdir,"03-Candidates-fl_dom.tsv")
hmmdom<-read_hmm_dom(hmmdomfile)
id.file<-file.path(outdir,"02-Candidates-all.id")
id<-fread(id.file)

#Filter hmm hits by evalue
hmm.sele<-unique(hmmdom[Evalue.Full<0.01,])

hmm.wide<-merge(dcast(hmm.sele,Target_name~Query_name,value.var = "Evalue.Full",fun.aggregate = length),id,by.x = "Target_name",by.y="Num",all.x = T)
hmm.wide[,Comp:=ifelse(HisKA>0,yes = 1,no = 0) +
                      ifelse(HATPase_c>0,yes = 1,no = 0) +
                      ifelse(Response_reg >0,yes = 1,no = 0)
                      ,by=.(Target_name)]


hmm.final<-hmm.wide[Comp>=2,]

#Check if there are any sequences that could not be found
list.missing<-setdiff(hmm.final$Target_name,fasta$ID)
if (length(list.missing)>0) {
  message("Some hmm hits could not be found in the fasta!") 
  message(paste(list.missing,sep = "\n"))
} 

#Select sequences in the fasta
fasta.sele<-fasta[ID %in% hmm.final$Target_name,]

write(fasta.sele$Text,file.path(outdir,"03-Candidates-filtered.fasta"))
