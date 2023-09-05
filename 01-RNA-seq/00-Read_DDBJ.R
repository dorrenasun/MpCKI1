rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update = T)

check_package<-function(pkg){
  if (!require(pkg,character.only = T,quietly = T)){
    # message(paste(pkg,"is not here. Install."))
    BiocManager::install(pkg,ask=F)
    library(pkg,character.only = T)
  }
}

pkg_list<-c("data.table","R.utils","rvest","XML")
for (p in pkg_list) {check_package(p)}

out.path<-"./Salmon_fem_Aki"
if (!dir.exists(out.path)) dir.create(out.path)

#read folders from sever home
homepage = "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA014/DRA014460/"
home<-read_html(homepage)
home.list<-html_text2(html_elements(home,"a"))
home.sub<-grep("/",home.list,value = T)

#Download xml tables
home.xml<-grep("xml",home.list,value = T)
for (x in home.xml){
  xml.url<-paste0(homepage,x)
  xml.file<-file.path(out.path,x)
  if (!file.exists(xml.file)) download.file(url=xml.url,destfile = xml.file)
}

#Download fasta files from each directory
# for (d in home.sub){
#   dname<-sub("\\/","",d)
#   out.path.sub<-file.path(out.path,dname)
#   if (!dir.exists(out.path.sub)) dir.create(out.path.sub)
#   #Find url for each file
#   sub.url<-paste0(homepage,d)
#   sub.text<-html_text2(html_elements(read_html(sub.url),"a"))
#   sub.files<-grep("fastq",sub.text,value = T)
#   for (f in sub.files){
#     url.file<-paste0(sub.url,f)
#     dest.file<-file.path(out.path.sub,f)
#     if (!file.exists(dest.file)) download.file(url=url.file,destfile = dest.file)
#   }
# }

#Process message from xml tables
file.run<-list.files(out.path,pattern = "run.xml",full.names = T )
run.top <- xmlRoot(xmlParse(file.run)) #root node
dt.run<-data.table()
for (i in 1:xmlSize(run.top)) {                
  alias<-xmlGetAttr(node = run.top[[i]],    
                       name = "alias")
  title<-xmlValue( run.top[[i]][["TITLE"]])
  center<-xmlGetAttr(node = run.top[[i]],    
                     name = "center_name")
  acc<-xmlGetAttr(node = run.top[[i]][["EXPERIMENT_REF"]],    
                  name = "accession")
  ent.run<-as.data.table(as.list(c(alias=alias,
                                   accession=acc,
                                   center_name=center,
                                   title=title
                                   )))
  dt.run<-rbind(dt.run,ent.run)
}
dt.run[,PRIMARY_ID:=sub(".*(SAM.*)","\\1",title)]
file.sample<-list.files(out.path,pattern = "sample.xml",full.names = T )
sample.top <- xmlRoot(xmlParse(file.sample)) #root node
dt.sample<-data.table()
for (i in 1:xmlSize(sample.top)) {
  id<-xmlToList(sample.top[[i]][["IDENTIFIERS"]],addAttributes = F)
  name<-xmlToList(sample.top[[i]][["SAMPLE_NAME"]])
  attributes <- xmlToDataFrame(sample.top[[i]][["SAMPLE_ATTRIBUTES"]])
  attributes.ls<-setNames(as.list(attributes$VALUE),attributes$TAG)
  
  ent.sample<-as.data.table(c(id,name,attributes.ls))
  dt.sample<-rbind(dt.sample,ent.sample)
}

dt.info<-merge(dt.run,dt.sample,by="PRIMARY_ID",all=T)
write.table(dt.info,file=file.path(out.path,"RunInfo.tsv"),sep="\t",col.names=T,row.names=F,quote=F)
