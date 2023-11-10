#################################### Environment ###############################
rm(list = ls())
message("Running 01-RetrieveFasta.R !")

if (!require(readxl)) install.packages("readxl")
if (!"readxl" %in% .packages()) library(readxl)

Round=1
source("00-Functions.R")

#Copy scripts
thisfile<-file.path(codedir_ori,"01-RetrieveFasta.R")
myCopy(thisfile,codedir)

#################################### Parameter ###############################

cutoff=50

#################################### Main ###############################
#Retrieve blast sequences from annotated genomes
{
  message("\nRetrieving fasta for genome hits!")
  
  #Copy input files
  inFile.tsv.Genome<-file.path(outdir_ori,"01-AnnotatedGenomes-blast.tsv")
  inFile.list.Genome<-file.path(indir_ori,"List_of_Genomes.xlsx")
  myCopy(inFile.tsv.Genome,outdir)
  myCopy(inFile.list.Genome,indir)
  
  #Large input files - do not copy
  inFile.fasta.Genome<-list.files(path = Genomedir,pattern = "merged.fasta",full.names = T)
  inFile.ID2Src.Genome<-list.files(path = Genomedir,pattern = "ID2Source",full.names = T)
  
  #Load blast results
  tsv.Genome<-fread(inFile.tsv.Genome,header = F)
  names(tsv.Genome)<-unlist(strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos",split=" "))

  #Correct "SM_"issues caused by blast
  list.SM<-unique(tsv.Genome[grep(saccver,pattern = "SM_",fixed = T),saccver])
  if (length(list.SM)>0){
    message(paste0("Correct names for:\n",paste0(list.SM,collapse = "\n" )))
    tsv.Genome[saccver %in% list.SM,saccver:=sub("SM_","Sm_",saccver)]
  }
  
  #Correct "ADC"issues caused by blast
  list.ADC<-unique(tsv.Genome[grep(saccver,pattern = "ADC",fixed = T),saccver])
  if (length(list.ADC)>0){
    message(paste0("Correct names for:\n",paste0(list.ADC,collapse = "\n" )))
    tsv.Genome[saccver %in% list.ADC,saccver:=sub("ADC","Adc",saccver)]
  }
  
  #Load ID2Src
  message("Loading ID2Src file for genomes...")
  ID2Src.Genome<-fread(inFile.ID2Src.Genome)
  #merge with ID2Src
  tsv.Genome<-merge(tsv.Genome,ID2Src.Genome,by.x = "saccver",by.y = "pID",all.x = T)
  #Check any IDs not found in ID2Src
  if (sum(is.na(tsv.Genome$Source))>0){
    message("Some blast hits not found in Genome ID2Src! Problem with parsing")
  }
   
  #Load list of Genomes to include  
  list.Genome<-as.data.table(read_excel(inFile.list.Genome,sheet = 1 ))
  #Check if there are any unmatches in Source names
  list.unmatch<-setdiff(ID2Src.Genome$Source,list.Genome$`File name`)
  if (length(list.unmatch)>0) {
    message("Unmatching sources! Check your genome parsing!")
    message(paste(list.unmatch,sep = "\n"))
  }  
  list.Src<-list.Genome[(!is.na(Include))&(Include=="T"),`File name`]

  #Filter ID by selected genome
  tsv.Genome.sele<-tsv.Genome[Source %in% list.Src]

  #Select the top N genes from each species
  setkeyv(tsv.Genome.sele,c("qaccver","Source","evalue"))
  tsv.Genome.sele[,No:=match(Gene,unique(Gene)),by=c("qaccver","Source")]
  tsv.Genome.cut<-tsv.Genome.sele[No<=cutoff,]
  
  #Keep only the best hit from the primary transcript for each gene
  setkeyv(tsv.Genome.cut,c("Gene","var","evalue"))
  tsv.Genome.cut[,primary:=(var==min(var) & !duplicated(var)),by=.(Gene)]
  tsv.Genome.final<-tsv.Genome.cut[primary==T,]
  
  #Load fasta sequences
  message("Loading fasta sequence for genomes...")
  fasta.Genome<-read_fasta(inFile.fasta.Genome)
  #Check if there are any unmatches in blast hits
  list.missing<-setdiff(tsv.Genome.final$saccver, fasta.Genome$ID)
  if (length(list.missing)>0) {
    message("Unmatching IDs! Some Genome blast hits not found!") 
    message(paste(list.missing,sep = "\n"))
  } 
  
  fasta.Genome.sele<-fasta.Genome[ID %in% tsv.Genome.final$saccver,]
 
  write(fasta.Genome.sele$Text,file.path(outdir,"01-AnnotatedGenomes-blast_hits.fasta"))
  message("Fasta retrieved for genome hits!")
  
}

# Retrieve OneKP results
{
  message("\nRetrieving fasta for onekp hits!")
  
  #Copy input files
  inFile.tsv.OneKP<-file.path(outdir_ori,"01-OneKP-Nonseed-blast.tsv")
  inFile.list.OneKP<-file.path(indir_ori,"List_of_OneKP.xlsx")
  myCopy(inFile.tsv.OneKP,outdir)
  myCopy(inFile.list.OneKP,indir)
  
  #Large input files - do not copy
  inFile.fasta.OneKP<-list.files(path = OneKPdir,pattern = "merged.fasta",full.names = T)

  #Load blast results
  tsv.OneKP<-fread(inFile.tsv.OneKP,header = F)
  
  names(tsv.OneKP)<-unlist(strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos",split=" "))
 
  
  #Load list of OneKP items to include 
  list.OneKP<-as.data.table(read_excel(inFile.list.OneKP,sheet = 1 ))
  #Check if there are any unmatches in OneKP codes
  tsv.OneKP[,Code:=sub(".*_([A-Z]+)-[0-9].*","\\1",saccver)]
  list.unmatch<-setdiff(tsv.OneKP$Code, list.OneKP$Code)
  if (length(list.unmatch)>0) warning("Unmatching Codes! Check your OneKP parsing!")  
  
  
  #Filter by OneKP sources and ranking in the same species
  list.OneKP.sele<-list.OneKP[(!is.na(Include))&(Include=="T"),Code]
  tsv.OneKP.sele<-tsv.OneKP[Code %in% list.OneKP.sele,]
  setkeyv(tsv.OneKP.sele,c("qaccver","Code","evalue"))
  tsv.OneKP.sele[,No:=seq(1:.N),by=c("qaccver","Code")]
  tsv.OneKP.cut<-tsv.OneKP.sele[No<=cutoff,]
  
  #Keep the best hit for each gene
  setkeyv(tsv.OneKP.cut,c("saccver","evalue"))
  tsv.OneKP.final<-tsv.OneKP.cut[!duplicated(saccver)]
  
  #Load OneKP fasta
  message("Loading fasta file for OneKP...")
  fasta.OneKP<-read_fasta(inFile.fasta.OneKP)
  #Check if there are any missing IDs
  list.missing<-setdiff(tsv.OneKP.final$saccver, fasta.OneKP$ID)
  if (length(list.missing)>0) {
    message("Unmatching IDs! Some OneKP blast hits not found!")  
    message(paste0(list.missing,sep = "\n"))
    }
  
  fasta.OneKP.sele<-fasta.OneKP[ID %in% tsv.OneKP.final$saccver,]
  write(fasta.OneKP.sele$Text,file.path(outdir,"01-OneKP-blast_hits.fasta"))
  message("Fasta retrieved for onekp hits!")
  
    
}

# Retrieve Dong2022 results
{
  message("\nRetrieving fasta for Dong2022 hits!")
  
  #Copy input files
  inFile.tsv.Dong<-file.path(outdir_ori,"01-Dong-blast.tsv")
  inFile.list.Dong<-file.path(indir_ori,"List_of_Dong2022.xlsx")
  myCopy(inFile.tsv.Dong,outdir)
  myCopy(inFile.list.Dong,indir)
  
  #Large input files - do not copy
  inFile.fasta.Dong<-list.files(path = Dongdir,pattern = "merged.fasta",full.names = T)
  inFile.ID2Src.Dong<-list.files(path = Dongdir,pattern = "ID2Source.txt",full.names = T)
  
  #Load blast results
  tsv.Dong<-fread(inFile.tsv.Dong,header = F)
  names(tsv.Dong)<-unlist(strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos",split=" "))
  
  #Load ID2Src
  message("Loading ID2Src file for Dong 2022...")
  ID2Src.Dong<-fread(inFile.ID2Src.Dong)
  #merge with ID2Src
  tsv.Dong<-merge(tsv.Dong,ID2Src.Dong,by.x = "saccver",by.y = "pID",all.x = T)
  #Check any IDs not found in ID2Src
  if (sum(is.na(tsv.Dong$Source))>0){
    message("Some blast hits not found in Dong 2022 ID2Src! Problem with parsing")
  }
  
  #Load list of Genomes to include  
  list.Dong<-as.data.table(read_excel(inFile.list.Dong,sheet = 1,range = cell_cols("A:I") ))
  #Check if there are any unmatches in Source names
  list.unmatch.Dong<-setdiff(ID2Src.Dong$Source,list.Dong$`File Name`)
  if (length(list.unmatch.Dong)>0) {
    message("Unmatching sources for Dong 2022! Check your genome parsing!")
    message(paste(list.unmatch.Dong,sep = "\n"))
  }  
  list.Src.Dong<-list.Dong[(!is.na(Include))&(Include=="T"),`File Name`]
  
  #Filter ID by selected transcriptome
  tsv.Dong.sele<-tsv.Dong[Source %in% list.Src.Dong,]
  

  #Recongize gene, isoforms and different proteins from TRINITY output
  tsv.Dong.sele[grepl(".*c[0-9]+_g[0-9]+_.*",saccver),`:=`(Gene=sub("_i.*\\.p.*","",saccver),
                                                           var=sub(".*_(i.*\\.p.*)","\\1",saccver)
  )]
  tsv.Dong.sele[is.na(var),`:=`(Gene=saccver,var=1)]
  
  #Select the top N genes from each species
  setkeyv(tsv.Dong.sele,c("qaccver","Source","evalue"))
  tsv.Dong.sele[,No:=match(Gene,unique(Gene)),by=.(qaccver,Source)]
  tsv.Dong.cut<-tsv.Dong.sele[No<=cutoff,]
  
  
  #Keep only the best hit from each gene
  setkeyv(tsv.Dong.cut,c("Gene","evalue"))
  tsv.Dong.cut[,primary:=(evalue==min(evalue)),by=.(qaccver,Gene)]
  tsv.Dong.final<-tsv.Dong.cut[primary==T,]
  
  
  #Load Dong fasta
  message("Loading fasta file for Dong 2022...")
  fasta.Dong<-read_fasta(inFile.fasta.Dong)
  #Check if there are any missing IDs
  list.missing<-setdiff(tsv.Dong.final$saccver, fasta.Dong$ID)
  if (length(list.missing)>0) {
    message("Unmatching IDs! Some Dong blast hits not found!")  
    message(paste(list.missing,sep = "\n"))
  }
  
  
  fasta.Dong.sele<-fasta.Dong[ID %in% tsv.Dong.final$saccver,]
  write(fasta.Dong.sele$Text,file.path(outdir,"01-Dong-blast_hits.fasta"))
  message("Fasta retrieved for Dong hits!")
  
  
}
