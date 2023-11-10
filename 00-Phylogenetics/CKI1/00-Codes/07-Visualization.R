#################################### Environment ###############################
rm(list=ls())
message("Running 07-Visualization.R !")
Round=2
source("00-Functions.R")

check_package<-function(pkg){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",update = T)
  if (!require(pkg,character.only = T,quietly = T)){
    # message(paste(pkg,"is not here. Install."))
    BiocManager::install(pkg,ask=F)
    library(pkg,character.only = T)
  }
}
pkg_list<-c("data.table","ggplot2","ggtree","phytools","treeio","ggmsa","dplyr","tidytree","ggnewscale","aplot","ggtext")
for (p in pkg_list) {check_package(p)}

#Copy scripts
thisfile<-file.path(codedir_ori,"06-ExtractFromSubtree.R")
myCopy(thisfile,codedir,overwrite = T)
############################# Plot Settings ####################################
#Parameters
lwd <- 1/.pt
fontsize <- 6
colors_Clade<-c(
  "tan4","#8DD35F","#2CA02C","#88AA00","#AB37C8","#5F5FD3","#5599FF","#0044AA"
)

# order.clade<-c([1]"Chlorophytes",
#                [2]"Streptophyte algae",
#                [3]"Mosses",
#                [4]"Liverworts",
#                [5]"Hornworts",
#                [6]"Lycophytes",
#                [7]"Ferns",
#                [8]"Gymnosperms",
#                [9]"Angiosperms",
#                [10]"Basal Angiosperms",
#                [11]"Monocots",
#                [12]"Magnoliids",
#                [13]"Eudicots",
#                [14]"Unknown")

names(colors_Clade)<-order.clade[c(2:9)]

colors_CKI1 <- c("CKI1" = "grey20","Others" = "grey50")

clustal.c1 <- "orange"
clustal.c2 <- "red3"
clustal.c3 <- "royalblue"
clustal.c4 <- "springgreen4"
        
fill.clustal <- c(
  "G" = clustal.c1,
  "P" = clustal.c1,
  "S" = clustal.c1,
  "T" = clustal.c1,
  "H" = clustal.c2,
  "K" = clustal.c2,
  "R" = clustal.c2,
  "F" = clustal.c3,
  "W" = clustal.c3,
  "Y" = clustal.c3,
  "I" = clustal.c4,
  "L" = clustal.c4,
  "M" = clustal.c4,
  "V" = clustal.c4
)
      
fill.domain <- c(
  # "black", #TM
  "#4096e5", #TM
  # "#0072B2", #TM
  # "#56B4E9", #TM
  "#ED7D31", #CHASE
  #"#CC79A7", #CHASE
  "#595959", #HK
  #"orange", #HK
  "#a6a6a6" #RR
  #"#009E73" #RR
)


#################################### Main ######################################
#Inputs & outputs
file.tree<-list.files(path = koutdir_ori,pattern = "08.*ordered\\.newick",full.names = T)
if (!file.exists(file.tree)){
  stop("Tree file not found in 03-KeyOutput! Exit!")
}else{
  myCopy(file.tree,koutdir,overwrite = T)
}

file.align<-list.files(path=outdir,"07-Subtree_[LE]INSI.fasta",full.names = T)
if (!file.exists(file.align)){
  stop("Alignment file not found in 02-Output! Exit!")
}

file.ipr<-file.path(outdir,"07-Subtree-IPR.tsv")
if (!file.exists(file.ipr)){
  source("07-ProcessIPR.R")
}
if (file.exists(file.ipr)) myCopy(file.ipr,koutdir)

file.names<-file.path(indir_ori,"Name_substitute.csv")
if (!file.exists(file.names)){
  message("File for alternative names not found in 01-Input.")
}else{
  myCopy(file.names,indir,overwrite = T)
}

file.id<-file.path(outdir,"06-Subtree.id")
if (!file.exists(file.id)){
  stop("File for sequence ids not found in 02-Output!Exit!")
}
    
#Load the ordered tree
tree1<-read.newick(file.tree)
tree1$tip.label<-gsub("[ ']","",tree1$tip.label)
tree1$node.label<-gsub("[ ']","",tree1$node.label)

#Extract tip names
annot<-data.table(No=1:length(tree1$tip.label),Name=tree1$tip.label)


#Read alternative names & species info from file
if (file.exists(file.names)){
  names<-fread(file.names)
  annot_sub<-merge(annot,names,by.x = "Name",by.y="Long_Name",all.x=T)
  annot_sub[is.na(Short_Name),Short_Name:=Name]
}else{
  annot_sub<-copy(annot)
  annot_sub[,Short_Name:=Name]
}


#Find IDs for target sequences
id<-fread(file.id)
id.sub<-id[,.(Name_sep=unlist(strsplit(Name,"[;\\|]")),Species,Clade,Tag),by=ID]

annot_sep<-annot_sub[,.(Name,Short_Name,Name_sep=unlist(strsplit(Name,"[;\\|]"))),by=No]
annot_sep1<-merge(annot_sep,id.sub,by="Name_sep",all.x = T)
annot_uniq<-unique(annot_sep1[,.(No,Name,Short_Name,label=ID,Species,Clade,Tag)])

#Check if all IDs are found
list.missing<-annot_uniq[is.na(label),Name]
if (length(list.missing)>0){
  message("Some IDs not found for sequence names!")
  message(paste0(list.missing,collapse = "\n"))
}

#Create label with species:
annot_uniq[,Full_Name:=paste0("[<i>",Species,"</i>] ",Name)]

#Fix order of legend:
annot_uniq[,Clade:=factor(Clade,levels = order.clade)]

#An empty column for generating dotted lines
annot_uniq[,Blank:=""]


#Substitue tree labels
tree1$tip.label<-annot_uniq[match(tree1$tip.label,Name),label]
#Convert to tibble
tree1.tbl<-as_tibble(tree1)

#Join info
tree2<-full_join(tree1,annot_uniq,by="label")
#Convert to tibble
tree2.tbl<-as_tibble(tree2)


{ #Extract various nodes
  node.SmCKI1L<-tree2.tbl$node[!is.na(tree2.tbl$Name)&grepl("Sm_111212",tree2.tbl$Name)]
  node.hwCKI1L<-tree2.tbl$node[!is.na(tree2.tbl$Name)&tree2.tbl$Name=="Apun_utg000041l.200.1"]
  node.CKI1L<-MRCA(tree1.tbl,node.SmCKI1L,node.hwCKI1L)$node
  
  node.AtCKI1<-tree2.tbl$node[!is.na(tree2.tbl$Name)& grepl("AT2G47430",tree2.tbl$Name)]
  node.MpCKI1<-tree2.tbl$node[!is.na(tree2.tbl$Name)& grepl("Mp2g02750",tree2.tbl$Name)]
  node.CKI1.strict<-MRCA(tree1.tbl,node.AtCKI1,node.MpCKI1)$node
  
  node.CKI1.all<-MRCA(tree1.tbl,node.SmCKI1L,node.AtCKI1)$node
  node.CKI1.parent<-parent(tree1.tbl,node.CKI1.all)$node
  sub.CKI1<-offspring(tree1.tbl, node.CKI1.all)
  
  node.SmAHK1<-tree2.tbl$node[!is.na(tree2.tbl$Name)&grepl("Sm_104037",tree2.tbl$Name)]
  node.AtAHK1<-tree2.tbl$node[!is.na(tree2.tbl$Name)&grepl("AT2G17820",tree2.tbl$Name)]
  node.AHK<-MRCA(tree1.tbl,node.SmAHK1,node.AtAHK1)$node
  
  node.AtAHK2<-tree2.tbl$node[!is.na(tree2.tbl$Name)&grepl("AT5G35750",tree2.tbl$Name)]
  node.MpCHK1<-tree2.tbl$node[!is.na(tree2.tbl$Name)&grepl("Mp2g03050",tree2.tbl$Name)]
  node.CHK<-MRCA(tree1.tbl,node.MpCHK1,node.AtAHK2)$node
  
}

#Add Core info to tree2.tbl
tree2.tbl$Tag[is.na(tree2.tbl$Name)&tree2.tbl$node %in% c(sub.CKI1$node,node.CKI1.all)]<-"Core"


#Read alignment file:
align<-read_fasta(file.align)
target_seq<-"FANASHDIRGALAG" # based on AtCKI1 seq
AtCKI1_seq<-align[ID==id.sub[grepl("AT2G47430",Name_sep),ID],Sequence]
target_start<-unlist(gregexpr(target_seq, AtCKI1_seq))
target_end<-target_start+nchar(target_seq)-1
msa<-as.data.table(tidy_msa(file.align,start = target_start,end=target_end))

# msa<-merge(annot_uniq,msa_ori,by.x="Acc",by.y = "name",all = T)

#Load IPR information
ipr<-fread(file.ipr,fill=T)
annot_len<-merge(annot_uniq,unique(ipr[,.(query,len,pos=0,wd=len/max(len))]),by.x="label",by.y = "query",all.x = T)
annot_ipr<-merge(annot_uniq,unique(ipr[,.(query,name,len,pos=((stop+start-1-len)/2)/max(len),wd=(stop-start+1)/max(len))]),by.x="label",by.y = "query",all.x = T)
annot_ipr[,name:=factor(name,levels = unique(name))]


{ #Plot the tree: with all info
  p1 <- ggtree(as.treedata(tree2.tbl), ladderize = F, color="black",size=0.5*lwd) +
    # ggplot(tree2,ladderize=F)+
    # geom_tree(aes(color=IsCKI,size=IsCKI),size=1*lwd)+
    geom_nodelab(
      mapping = aes(label = label),
      color="grey30",
      size = 3 * lwd,
      hjust = -0.2
    ) +
    scale_color_manual(values = colors_CKI1,guide = 'none') +
    # scale_size_manual(values = c("CKI1" = lwd * 1.5, "Others" = lwd * 1)) +
    ggnewscale::new_scale_colour()  +
    geom_tippoint(
      mapping = aes(x, y, color = Clade),
      size = 1.5,
      alpha = 1,
      shape = 15,
      position = position_nudge(x = 0.05, y = 0),
      inherit.aes = T
    ) +
    scale_color_manual(values = colors_Clade) +
    ggnewscale::new_scale_color()  +
    geom_tiplab(
      mapping = aes(color = Clade, label = Blank),
      align = T,
      offset = 0.5,
      size = fontsize * lwd,
      linetype = "dotted",
      linesize = 0.5 * lwd
    ) +
    scale_color_manual(values = colors_Clade,guide="none") +
    # theme_classic()+
    # theme_tree()+
    geom_treescale(
      x = 0,
      y = -max(annot_sub$No) + 1,
      fontsize = fontsize * lwd,
      linesize = 1 * lwd,
      offset = 1
    ) +

    geom_facet(
      panel = "Domain",
      data = annot_len[,.(label,len,wd)],
      geom = geom_tile,
      mapping = aes(x = 0,width=wd,height=0.7),
      fill = "white",#"grey30"#,
      color="grey50",
      linewidth = 0.3 * lwd
      ) +
    ggnewscale::new_scale_fill()  +
    geom_facet(
      panel = "Domain",
      # data = annot_dom_all[,.(label,Domain=Query_name,pos,wd)],
      data = annot_ipr[!name=="PRibTrfase_dom",.(label,Domain=name,pos,wd)],
      geom = geom_tile,
      mapping = aes(x = pos,width=wd,height=0.7,fill=Domain),
      color="grey50",
      linewidth = 0.3 * lwd
      # fill = "grey30"#,
        # linewidth = 2 * lwd
    ) +
    # scale_fill_manual(values = fill.domain)+
    ggnewscale::new_scale_fill()  +
    geom_facet(
      panel = "MSA",
      data = msa[, .(id = name, position, character)],
      geom = geom_tile,
      mapping = aes(x = position, fill = character),
      alpha = 0.6
    ) +
    scale_fill_manual(values = fill.clustal, na.value = "white",guide = 'none') +
    geom_facet(
      panel = "MSA",
      data = msa[, .(id = name, position, character)],
      geom = geom_text,
      mapping = aes(x = position, label = character),
      color = "grey30",
      size = 6 * lwd,
      vjust = "middle"
    ) +
    ggnewscale::new_scale_color()  +
    geom_facet(
      panel = "ID",
      data = annot_uniq[, .(id = label,pos = -0.05, clade=Clade)],
      geom = geom_point,
      mapping = aes(x = pos, color = clade),
      size = 1.5,
      alpha = 1,
      shape = 15
      # color = "grey30",
      # size = 6 * lwd,
      # vjust = "middle"
    ) +
  scale_color_manual(values = colors_Clade,guide = 'none')+
    ggnewscale::new_scale_color()+
    geom_facet(
      panel = "ID",
      data = rbind(annot_uniq[, .(id = label,pos = 0, tx=Full_Name,cl=Clade,cki=Tag)],
                   data.table(id = annot_uniq$label[1],pos = 1, tx="",cl=NA,cki="Others")),
      geom = geom_richtext,
      # parse = TRUE,
      mapping = aes(x=pos,label = tx,color=cki),
      label.color = NA,
      label.padding=unit(c(0, 0, 0, 0), "lines"),
      color = "black",
      # color="grey30",
      fill=NA,
      size = 6 * lwd,
      # fontface="italic",
      vjust = "middle",
      hjust="left"
    ) +
    # scale_color_manual(values = colors_CKI1,guide = 'none') +
    #Y axis labels
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(
      expand = expansion(0, 0.6),
      name = "IDs",
      breaks = unique(msa$No),
      labels = unique(msa$label),
      position = "right"
    ) +

    # theme_grey() +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.title = element_blank(),
      axis.text.y=element_blank(),
      # axis.text.y = element_text(
      #   family = "sans",
      #   size = fontsize,
      #   color = "black",
      #   face = "plain"
      # ),
      legend.text = element_text(size=fontsize,color = "black"),
      legend.title = element_text(size=fontsize,color = "black"),
      legend.key.size = unit(0.2,"cm"),
      legend.position = "left"
    )

  # facet_grid(~panel,space="free_x")
  
  p1 <- facet_widths(p1, c(1,0.6, 0.3,1))
plot(p1)
outfile1<-file.path(outdir,sub(".newick",".pdf",basename(file.tree),fixed = T))
ggsave(outfile1,p1,width = 22, height = 26,units = "cm")

}

stop()

{p2<-ggtree(tree2,ladderize=F, layout="circular",aes(color=Clade),size=0.5*lwd,open.angle = 178) +
    geom_tippoint(
      mapping = aes(x, y,angle, color = Clade),
      size = 0.8*lwd,
      alpha = 1,
      shape = 16,
      position = position_nudge(x = 0.04, y = 0),
      inherit.aes = T
    ) +
    geom_treescale(
      x = 1,
      y=0,
      color="grey30",
      # y = -max(annot_sub$No) + 1,
      fontsize = fontsize * lwd,
      linesize = 0.5 * lwd,
      offset = 1
    ) +
    geom_point2(aes(subset=(node==node.MpCKI1)),
                position = position_nudge(x = 0.08, y = 0),
                shape=16, size=0.8*lwd, color='red')+
    # geom_tiplab()+
    # geom_tiplab2(aes(subset=(node==node.MpCKI1),hjust=0,label="MpCKI1"))+
    scale_y_reverse()+
    geom_hilight(mapping=aes(subset = node %in% node.CKI1.strict),
                 fill="#d7eef4ff",
                 to.bottom = T,
                 extend=0.7,
                 # type = "gradient", gradient.direction = 'tr',
                 alpha = .8) +
    geom_cladelab(node=node.CKI1.strict, label="CKI1", angle=0, 
                  barcolor="#216778",
                  fontsize=8*lwd, offset=.7, vjust=.5)  + 
    geom_cladelab(node=node.CKI1L, label="CKI1L", angle=0, 
                  barcolor="#668000",
                  fontsize=8*lwd, offset=.2, vjust=.5)  + 
    geom_cladelab(node=node.AHK, label="AHK1", angle=0, 
                  barcolor="#786721",
                  fontsize=8*lwd, offset=.2, vjust=.5)  + 
    geom_cladelab(node=node.CKR, label="CHK", angle=0, 
                  barcolor="#AA4400",
                  fontsize=8*lwd, offset=.2, vjust=.5)  + 
    geom_hilight(mapping=aes(subset = node %in% node.CKI1L),
                 fill="#e0eccdff",
                 to.bottom = T,
                 extend=0.2,
                 # type = "gradient", gradient.direction = 'tr',
                 alpha = .8) +
    geom_hilight(mapping=aes(subset = node %in% node.AHK),
                 fill="#f4eed7ff",
                 to.bottom = T,
                 extend=0.2,
                 # type = "gradient", gradient.direction = 'tr',
                 alpha = .8) +
    geom_hilight(mapping=aes(subset = node %in% node.CKR),
                 fill="#faecd9ff",
                 to.bottom = T,
                 extend=0.2,
                 # type = "gradient", gradient.direction = 'tr',
                 alpha = .8) +
    
    # coord_polar(theta = "y",direction = 1,start=-pi/2,clip = "on")+
    # ggplot(tree2,ladderize=F)+
    # geom_tree(aes(color=IsCKI,size=IsCKI),size=1*lwd)+
    # coord_flip()+
    scale_color_manual(values = colors_Clade) +
    # theme_grey()#+
    
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          legend.key.size = unit(0.2,"cm"),
          legend.text = element_text(size=6,color = "grey30"),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.position = "right"

          )
  plot(p2)
  outfile1<-file.path(OutDir,"07-Subtree_EINSI_80_np-fanplot.pdf")
  ggsave(outfile1,p2,width = 12, height = 8,units = "cm")
  
}

