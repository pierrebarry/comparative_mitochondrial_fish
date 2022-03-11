#--------------------------------|#
#     Pipeline for MitoZ          #
#--------------------------------/#

#.libPaths("/home/labosea1/R/x86_64-pc-linux-gnu-library/3.6")
#.libPaths("/usr/local/lib/R/site-library")
#.libPaths("/usr/lib/R/site-library")
#.libPaths("/usr/lib/R/library")

library(apex)
library(Biostrings)
library(stringi)
library(ggtree)
library(tidyverse)


SPECIES=list.dirs("/DATA/sdb1/Pierre/MitoZ/",full.names=F,recursive=F)

data_reference=data.frame(SPECIES=c("Aboye",
                                    "Afall",
                                    "Cgale",
                                    "Cjuli",
                                    "Dlabr",
                                    "Dpunt",
                                    "Eencr",
                                    "Gnige",
                                    "Hgutt",
                                    "Lbude",
                                    "Lmorm",
                                    "Mmerl",
                                    "Msurm",
                                    "Peryt",
                                    "Scabr",
                                    "Scant",
                                    "Scine",
                                    "Spilc",
                                    "Ssard",
                                    "Styph"),
                          REF=c("AboyeLi6",
                                "AfallLi5",
                                "CgaleLi5",
                                "CjuliLi1",
                                "DlabrLi3",
                                "DpuntFa4",
                                "EencrMu6",
                                "GnigeLi3",
                                "HguttFa3",
                                "LbudeMu1",
                                "LmormLi4",
                                "MmerlLi2",
                                "MsurmLi6",
                                "PerytMu5",
                                "ScabrFa2",
                                "ScantFa2",
                                "ScineMu5",
                                "SpilcGa3",
                                "SsardMu1",
                                "StyphFa2"
                          ))

SPECIES=c(
  #Aboye
  #"Afall",
  #
  #"Cjuli",
  #"Dlabr",
  #"Dpunt",
  #"Eencr"
  #"Gnige"
  #"Hgutt",
  #"Lbude",
  #"Lmorm",
  #"Mmerl",
  #"Msurm",
  #"Peryt",
  #"Scabr",
  #"Scant"
  #"Scine",
  #"Spilc",
  #"Ssard",
  #"Styph"
)


#Cgale

for (sp in SPECIES){
  
  SAMPLES=list.files(paste("/DATA/sdb1/Pierre/MitoZ/",sp,sep=""),full.names=F,recursive=F)
  
  SAMPLES=SAMPLES[SAMPLES!="AfallLi1" & 
                    SAMPLES!="HguttGa11" & SAMPLES!="HguttGa6" & SAMPLES!="HguttGa7" & SAMPLES!="HguttGa8" & SAMPLES!="HguttGa9" & SAMPLES!="HguttMu1" & SAMPLES!="HguttMu4" &
                    SAMPLES!="HguttGa4" &
                    SAMPLES!="LbudeMu5" & SAMPLES!="LbudeMu6"
                  & SAMPLES!="DpuntMu5" &
                    SAMPLES!="EencrMu4" &
                    SAMPLES!="ScantLi1" &
                    SAMPLES!="StyphFa3" & SAMPLES!="StyphFa4" & SAMPLES!="StyphGa1" & SAMPLES!="StyphGa2" &
                    SAMPLES!="StyphGa6" & SAMPLES!="StyphLi5" & SAMPLES!="StyphLi6" & SAMPLES!="StyphMu5"]
  
  # Extract reference
  system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref.fa",sep=""))
  
  setwd(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,sep=""))
  
  #if (sp=="Cgale"){
  #  
  #  system(paste("grep -w 'Ck141_593551' -A 1 ref.fa > ref_two.fasta",sep=""))
  #  system("mv ref_two.fasta ref.fa")
    
  if (sp=="Mmerl"){
    
    system(paste("grep -w 'Ck141_672172' -A 1 ref.fa > ref_two.fasta",sep=""))
    system("mv ref_two.fasta ref.fa")
    
  } else if (sp=="Scant") {
    
    system(paste("grep -w 'Ck141_600518' -A 1 ref.fa > ref_two.fasta",sep=""))
    system("mv ref_two.fasta ref.fa")    
    
  }
  
  for (sam in SAMPLES){
    
    # Reorder
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",sam,"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta tmp_subject.fasta",sep=""))
    system(paste("blastn -query ref.fa -subject tmp_subject.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    tmp<-read.table("tmp")
    if (sp=="Aboye"){
      tmp=tmp[tmp$V3>=80,]
      
    } else {
      tmp=tmp[tmp$V3>=90,]
    }
    
    tmp=tmp[order(tmp$V3,decreasing = T),]
    
    if (sam=="CgaleMu1"){
      system(paste("grep -w '","Ck141_447230","' -A 1 tmp_subject.fasta > tmp_subject_two.fasta",sep=""))
      
    } else if (sam=="CgaleMu5"){
      system(paste("grep -w '","Ck141_744009","' -A 1 tmp_subject.fasta > tmp_subject_two.fasta",sep=""))

    } else if (sam=="CgaleMu6"){
      system(paste("grep -w '","Ck141_320582","' -A 1 tmp_subject.fasta > tmp_subject_two.fasta",sep=""))
      
    } else if (sam=="CgaleMu2"){
      system(paste("grep -w '","Ck141_208874","' -A 1 tmp_subject.fasta > tmp_subject_two.fasta",sep=""))
      
    } else if (sam=="CgaleMu4"){
      system(paste("grep -w '","Ck141_316901","' -A 1 tmp_subject.fasta > tmp_subject_two.fasta",sep=""))
      
    } else {
      
      system(paste("grep -w '",strsplit(as.character(tmp$V2[1]),";")[[1]][1],"' -A 1 tmp_subject.fasta > tmp_subject_two.fasta",sep=""))
      
    }
    
    
    #if (sam=="CgaleMu1"){
    #  tmp=tmp[order(tmp$V4,decreasing = T),]
    #}
    
    if (tmp$V7[1]>tmp$V8[1]){
      
      if (sam=="CgaleFa2" | sam=="CgaleFa5"){
        tmp=tmp[order(tmp$V4,decreasing=T),]
        vect<-c()
        file<-read.delim("tmp_subject_two.fasta",header=F)
        vect<-c(vect,
                as.character(file[1,]),
                as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper"))
        )
        fileConn<-file("tmp_subject_two.fasta")
        writeLines(vect,fileConn)
        close(fileConn)
      } else if (sam=="CgaleMu1") {
        tmp=tmp[order(tmp$V5),]
        tmp=tmp[-2,]
        vect<-c()
        file<-read.delim("tmp_subject_two.fasta",header=F)
        vect<-c(vect,
                as.character(file[1,]),
                paste(as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper")),
                      as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[2,]$V8,tmp[2,]$V7)),case="upper")),
                      as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[3,]$V8,tmp[3,]$V7)),case="upper")),
                      collapse = NULL,sep="")
        )
        fileConn<-file("tmp_subject_two.fasta")
        writeLines(vect,fileConn)
        close(fileConn)
        
      } else if (sam=="CgaleGa4"){
        vect<-c()
        file<-read.delim("tmp_subject_two.fasta",header=F)
        vect<-c(vect,
                as.character(file[1,]),
                paste(as.character(reverseComplement(DNAString(substr(as.character(file[2,]),843,16487)),case="upper")),
                      collapse = NULL,sep="")
        )
        fileConn<-file("tmp_subject_two.fasta")
        writeLines(vect,fileConn)
        close(fileConn)
        
      } else {
        ## A SURVEILLER
        tmp=tmp[order(tmp$V5),]
        ## A SURVEILLER
        vect<-c()
        file<-read.delim("tmp_subject_two.fasta",header=F)
        if (nrow(tmp)==1){
          vect<-c(vect,
                  as.character(file[1,]),
                  as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper"))
          )
        } else if (nrow(tmp)==2) {
          vect<-c(vect,
                  as.character(file[1,]),
                  paste(as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper")),
                        as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[2,]$V8,tmp[2,]$V7)),case="upper")),
                        collapse = NULL,sep="")
          )
        } else {
          vect<-c(vect,
                  as.character(file[1,]),
                  paste(as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper")),
                        as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[2,]$V8,tmp[2,]$V7)),case="upper")),
                        as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[3,]$V8,tmp[3,]$V7)),case="upper")),
                        collapse = NULL,sep="")
          )        
        }
        
        fileConn<-file("tmp_subject_two.fasta")
        writeLines(vect,fileConn)
        close(fileConn)
      }

    } else {
      

        tmp=tmp[order(tmp$V7,decreasing = T),]
        

    }
    
    system("rm tmp_subject.fasta")
    system("mv tmp_subject_two.fasta tmp_subject.fasta")
    system(paste("python3 /home/labosea1/MitoZ/ANALYSIS/Mitogenome_reorder.py -f tmp_subject.fasta -r ref.fa",sep=""))
    file.rename("tmp_subject.fasta.reorder", paste(sam,"_mitoz.fasta",sep=""))
    
    system(paste("blastn -query ref.fa -subject ",sam,"_mitoz.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    tmp<-read.table("tmp")
    
    if (nrow(tmp)>1){
      
      if (grepl("Cgale",sam)){
        
        vect<-c()
        file<-read.delim(paste(sam,"_mitoz.fasta",sep=""),header=F)
        tmp=tmp[order(tmp$V4,decreasing=T),]
        vect<-c(vect,
                as.character(file[1,]),
                paste(substr(as.character(file[2,]),tmp[1,]$V7,tmp[1,]$V8),
                      collapse = NULL,sep=""))
        fileConn<-file(paste(sam,"_mitoz.fasta",sep=""))
        writeLines(vect,fileConn)
        close(fileConn)
      } else {
        
        vect<-c()
        file<-read.delim(paste(sam,"_mitoz.fasta",sep=""),header=F)
        tmp=tmp[order(tmp$V5),]
        
        if (sp=="Aboye"){
          vect<-c(vect,
                  as.character(file[1,]),
                  paste(substr(as.character(file[2,]),tmp[1,]$V7,tmp[1,]$V8),
                        substr(as.character(file[2,]),tmp[2,]$V7,tmp[2,]$V8),
                        substr(as.character(file[2,]),tmp[3,]$V7,tmp[3,]$V8),
                        collapse = NULL,sep=""))
        } else {
          
          vect<-c(vect,
                  as.character(file[1,]),
                  paste(substr(as.character(file[2,]),tmp[1,]$V7,tmp[1,]$V8),
                        substr(as.character(file[2,]),tmp[2,]$V7,tmp[2,]$V8),
                        collapse = NULL,sep=""))
        }

        
        

        
      }
      
      
    }
    
    system("rm tmp_subject.fasta")
    
    
  }
  
  #outrgroup
  run==0
  if (run==1){
    system(paste("blastn -query ref.fa -subject ../outgroup.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    tmp<-read.table("tmp")
    tmp=tmp[order(tmp$V5),]
    
    if (tmp$V7[1]>tmp$V8[1]){
      vect<-c()
      file<-read.delim("../outgroup.fasta",header=F)
      if (nrow(tmp)==1){
        vect<-c(vect,
                paste(as.character(file[1,]),"; circular",sep=" "),
                as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper"))
        )
      } else if (nrow(tmp)==2){
        vect<-c(vect,
                paste(as.character(file[1,]),"; circular",sep=" "),
                paste(as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper")),
                      as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[2,]$V8,tmp[2,]$V7)),case="upper")),
                      collapse = NULL,sep="")
        )
      } else {
        vect<-c(vect,
                paste(as.character(file[1,]),"; circular",sep=" "),
                paste(as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[1,]$V8,tmp[1,]$V7)),case="upper")),
                      as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[2,]$V8,tmp[2,]$V7)),case="upper")),
                      as.character(reverseComplement(DNAString(substr(as.character(file[2,]),tmp[3,]$V8,tmp[3,]$V7)),case="upper")),
                      collapse = NULL,sep="") 
        )
      }
      
      fileConn<-file("outgroup_reorder.fasta")
      writeLines(vect,fileConn)
      close(fileConn)
    } else {
      tmp=tmp[order(tmp$V7,decreasing = T),]
    }
    
    system(paste("python3 /home/labosea1/MitoZ/ANALYSIS/Mitogenome_reorder.py -f outgroup_reorder.fasta -r ref.fa",sep=""))
    file.rename("tmp_subject.fasta.reorder", paste(sam,"_mitoz.fasta",sep=""))
    system("rm tmp_subject.fasta")
    
    system(paste("blastn -query ref.fa -subject ",sam,"_mitoz.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    tmp<-read.table("tmp")
  }
  
  
  # Align 
  dir.create("Align")
  vect<-c()
  fasta_files=list.files()[grep(sp,list.files())]
  for (i in fasta_files){
    
    if (i=="CgaleMu6_mitoz.fasta"){
      
      vect<-c(vect,
              paste(">",strsplit(i,"_")[[1]][1],sep=""),
              as.character(file[2,]))
    } else {
      
      file<-read.delim(i,header=F)
      system(paste("blastn -query ref.fa -subject ",i," -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
      tmp<-read.table("tmp")
      if (tmp$V7>tmp$V8){
        vect<-c(vect,
                paste(">",strsplit(i,"_")[[1]][1],sep=""),
                as.character(file[2,])
                #as.character(reverseComplement(DNAString(stri_reverse(as.character(file[2,])))))
        )
      } else {
        vect<-c(vect,
                paste(">",strsplit(i,"_")[[1]][1],sep=""),
                as.character(file[2,]))
      }
    }
    

    
    
  }
  
  fileConn<-file(paste("Align/",sp,"_mitochondrial.fasta",sep=""))
  writeLines(vect,fileConn)
  close(fileConn)
  system(paste("/usr/local/anaconda/bin/mafft Align/",sp,"_mitochondrial.fasta > Align/",sp,"_mitochondrial_align.fa",sep=""))
  
  beeData<-read.multiFASTA(paste("Align/",sp,"_mitochondrial_align.fa",sep=""))
  pdf(paste("Align/",sp,"_align_data.pdf",sep=""),width=30,height=6)
  plot(beeData)
  dev.off()
  
  #fasttree
  #system(paste("/usr/local/anaconda/bin/fasttree -nt -gtr sequence_",gene[g],"_",name_species,"_align > tree_",gene[g],"_",name_species,sep=""))
  
  # iqtree
  dir.create("tree")
  system(paste("/usr/local/anaconda/bin/trimal -in Align/",sp,"_mitochondrial_align.fa -out tree/",sp,"_mitochondrial_align_trim.fa -htmlout tree/",sp,"_mitochondrial_align_trim.html -gt 1",sep=""))
  #system(paste("/usr/local/anaconda/bin/iqtree -s ",sp,"_mitochondrial_align_trim.fa -m MF -T 4",sep=""))
  system(paste("/usr/local/anaconda/bin/iqtree -s tree/",sp,"_mitochondrial_align_trim.fa -m TN+F+I -T 1 -redo",sep=""))
  system(paste("/usr/local/anaconda/bin/raxmlHPC -m GTRGAMMA -s tree/",sp,"_mitochondrial_align_trim.fa -n ",sp," -p 12345",sep=""))
  system(paste("mv RAxML_* tree/.",sep=""))
  
  tree <- read.tree(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/tree/",sp,"_mitochondrial_align_trim.fa.treefile",sep=""))
  image_info=data.frame(ind=tree$tip.label,
                        loc=substr(tree$tip.label,6,7),
                        sp=substr(tree$tip.label,1,5))
  
  image_info$loc=factor(image_info$loc,levels=c("Li","Mu","Fa","Ga"),
                        labels=c("Gulf of Lion",
                                 "Costa Calida",
                                 "Algarve",
                                 "Bay of Biscay")
  )
  p1<-ggtree(tree,col="black",size=1,alpha=0.5) %<+% image_info + 
    geom_tiplab(hjust=-0.15,size=3,fontface="bold",aes(color=loc)) +
    geom_tippoint(aes(color=loc),size=4,alpha=0.75) +
    theme_tree2() +
    scale_color_manual(name="Location",
                       values=RColorBrewer::brewer.pal(n=4,name="RdBu"))+
    scale_fill_manual(name="Species",values=wesanderson::wes_palette("Zissou1",n=length(levels(factor(image_info$sp))),type="continuous")) +
    #xlim(c(0,max(tree$edge.length))) +
    theme(legend.position = "bottom",
          plot.title=element_text(hjust=0.5)) +
    labs(title="iqtree")
  
  tree <- read.tree(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/tree/RAxML_bestTree.",sp,sep=""))
  image_info=data.frame(ind=tree$tip.label,
                        loc=substr(tree$tip.label,6,7),
                        sp=substr(tree$tip.label,1,5))
  
  image_info$loc=factor(image_info$loc,levels=c("Li","Mu","Fa","Ga"),
                        labels=c("Gulf of Lion",
                                 "Costa Calida",
                                 "Algarve",
                                 "Bay of Biscay")
  )
  p2<-ggtree(tree,col="black",size=1,alpha=0.5) %<+% image_info + 
    geom_tiplab(hjust=-0.15,size=3,fontface="bold",aes(color=loc)) +
    geom_tippoint(aes(color=loc),size=4,alpha=0.75) +
    theme_tree2() +
    scale_color_manual(name="Location",
                       values=RColorBrewer::brewer.pal(n=4,name="RdBu"))+
    scale_fill_manual(name="Species",values=wesanderson::wes_palette("Zissou1",n=length(levels(factor(image_info$sp))),type="continuous")) +
    #xlim(c(0,max(tree$edge.length)+max(tree$edge.length)*0.15)) +
    theme(legend.position = "bottom",
          plot.title=element_text(hjust=0.5)) +
    labs(title="RAxML")
  
  pdf(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/tree/",sp,"_tree.pdf",sep=""),width=30,height=7.5)
  print(ggpubr::ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend = "bottom"))
  dev.off()
  
  pdf(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/tree/",sp,"_tree_iqtree.pdf",sep=""),width=15,height=10)
  print(p1)
  dev.off()
  
  pdf(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/tree/",sp,"_tree_raxml.pdf",sep=""),width=15,height=10)
  print(p2)
  dev.off()
  
  
  # PopGen statistics
  
  system(paste("sh /home/labosea1/MitoZ/ANALYSIS/change_fasta.sh Align/",sp,"_mitochondrial_align.fa Align/",sp,"_mitochondrial_align_2.fa",sep=""))
  system(paste("mv Align/",sp,"_mitochondrial_align_2.fa Align/",sp,"_mitochondrial_align.fa",sep=""))
  
  dir.create("Pop_Gen")
  dir.create("VCF")
  system(paste("/usr/local/anaconda/bin/python3 /home/labosea1/MitoZ/ANALYSIS/run_mitoz.py ",sp," /home/labosea1/MitoZ/ANALYSIS/",sp,"/Align/",sp,"_mitochondrial_align.fa whole_mitochondrial",sep=""))
  system(paste("rm VCF/vcf_",sp,".vcf",sep=""))
  
  # PopGen statistics per gene
  #gene_file<-read.delim(paste("/DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/summary.txt",sep=""))
  #seq=strsplit(as.character(gene_file[1,])," ")[[1]][1]
  dir.create("GENE")
  
  list_gene=c()
  type_gene=c()
  start=c()
  stop=c()
  
  #if (sp=="Cgale"){
  #  
  #  system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_cds.fa",sep=""))
  #  system(paste("grep -w 'Ck141_593551' -A 1 ref_cds.fa > ref_cds_two.fasta",sep=""))
  #  system("mv ref_cds_two.fasta ref_cds.fa")
  #  system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_cds.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
  #  
  if (sp=="Mmerl"){
    
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_cds.fa",sep=""))
    system(paste("grep -w 'Ck141_672172' -A 1 ref_cds.fa > ref_cds_two.fasta",sep=""))
    system("mv ref_cds_two.fasta ref_cds.fa")
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_cds.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
    
  } else if (sp=="Scant") {
    
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_cds.fa",sep=""))
    system(paste("grep -w 'Ck141_600518' -A 1 ref_cds.fa > ref_cds_two.fasta",sep=""))
    system("mv ref_cds_two.fasta ref_cds.fa")
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_cds.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
    
  } else {
    
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
  }
  
  tmp<-read.table("tmp")
  tmp=tmp[tmp$V1==as.character(data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))]),]
  
  for (i in 1:nrow(tmp)){
    
    if (strsplit(as.character(tmp$V2[i]),";")[[1]][2] %in% list_gene == F){
      
      list_gene=c(list_gene,
                  strsplit(as.character(tmp$V2[i]),";")[[1]][2])
      type_gene=c(type_gene,
                  "CDS")
      
      if (tmp[i,]$V7<tmp[i,]$V8){
        start=c(start,
                as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][1])))
        stop=c(stop,
               as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][2])))
      } else {
        start=c(start,
                as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][2])))
        stop=c(stop,
               as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][1])))
      }
      
      
    }
    
    
  }
  
  #if (sp=="Cgale"){
  #  
  #  system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.rrna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_rrna.fa",sep=""))
  #  system(paste("grep -w 'Ck141_593551' -A 1 ref_rrna.fa > ref_rrna_two.fasta",sep=""))
  #  system("mv ref_rrna_two.fasta ref_rrna.fa")
  #  system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_rrna.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
  #  
  if (sp=="Mmerl"){
    
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.rrna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_rrna.fa",sep=""))
    system(paste("grep -w 'Ck141_672172' -A 1 ref_rrna.fa > ref_rrna_two.fasta",sep=""))
    system("mv ref_rrna_two.fasta ref_rrna.fa")
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_rrna.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
    
  } else if (sp=="Scant") {
    
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.rrna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_rrna.fa",sep=""))
    system(paste("grep -w 'Ck141_600518' -A 1 ref_rrna.fa > ref_rrna_two.fasta",sep=""))
    system("mv ref_rrna_two.fasta ref_rrna.fa")
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_rrna.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
    
  } else {
    
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.rrna -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
  }
  
  tmp<-read.table("tmp")
  tmp=tmp[tmp$V1==as.character(data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))]),]
  
  for (i in 1:nrow(tmp)){
    
    if (strsplit(as.character(tmp$V2[i]),";")[[1]][2] %in% list_gene == F){
      
      list_gene=c(list_gene,
                  strsplit(as.character(tmp$V2[i]),";")[[1]][2])
      type_gene=c(type_gene,
                  "rRNA")
      if (tmp[i,]$V7<tmp[i,]$V8){
        start=c(start,
                as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][1])))
        stop=c(stop,
               as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][2])))
      } else {
        start=c(start,
                as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][2])))
        stop=c(stop,
               as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][1])))
      }
    }
    
    
    
  }
  
  #if (sp=="Cgale"){
  #  
  #  system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.trna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_trna.fa",sep=""))
  #  system(paste("grep -w 'Ck141_593551' -A 1 ref_trna.fa > ref_trna_two.fasta",sep=""))
  #  system("mv ref_trna_two.fasta ref_trna.fa")
  #  system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_trna.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
  #  
  if (sp=="Mmerl"){
    
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.trna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_trna.fa",sep=""))
    system(paste("grep -w 'Ck141_672172' -A 1 ref_trna.fa > ref_trna_two.fasta",sep=""))
    system("mv ref_trna_two.fasta ref_trna.fa")
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_trna.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
    
  } else if (sp=="Scant") {
    
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.trna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_trna.fa",sep=""))
    system(paste("grep -w 'Ck141_600518' -A 1 ref_trna.fa > ref_trna_two.fasta",sep=""))
    system("mv ref_trna_two.fasta ref_trna.fa")
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject ref_trna.fa -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
    
  } else {
    
    system(paste("blastn -query Align/",sp,"_mitochondrial_align.fa -subject /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.trna -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
    
  }
  
  
  tmp<-read.table("tmp")
  tmp=tmp[tmp$V1==as.character(data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))]),]
  
  for (i in 1:nrow(tmp)){
    
    if (strsplit(as.character(tmp$V2[i]),";")[[1]][2] %in% list_gene == F){
      
      list_gene=c(list_gene,
                  strsplit(as.character(tmp$V2[i]),";")[[1]][2])
      type_gene=c(type_gene,
                  "tRNA")
      if (tmp[i,]$V7<tmp[i,]$V8){
        start=c(start,
                as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][1])))
        stop=c(stop,
               as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][2])))
      } else {
        start=c(start,
                as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][2])))
        stop=c(stop,
               as.numeric(gsub("\\D","",strsplit(strsplit(as.character(tmp$V2[i]),";")[[1]][4],":")[[1]][1])))
      }
      
    }
    
    
  }
  
  
  #system(paste("blastn -query /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta -subject /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
  #tmp<-read.table("tmp")
  #file<-read.delim(paste("/DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta",sep=""),header=F)
  #as.character(reverseComplement(DNAString(toupper(substr(as.character(file[2,]),tmp$V5[1],tmp$V6[1]))),case="upper"))
  
  for (g in list_gene){
    
    if (type_gene[which(list_gene==g)]=="CDS"){
      system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_cds.fa",sep=""))
      system(paste("grep -w '",g,"' -A 1 ref_cds.fa > ref_gene.fasta",sep=""))
    } else if (type_gene[which(list_gene==g)]=="rRNA"){
      system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.rrna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_rrna.fa",sep=""))
      system(paste("grep -w '",g,"' -A 1 ref_rrna.fa > ref_gene.fasta",sep=""))
    } else {
      system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.trna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_trna.fa",sep=""))
      system(paste("grep -w '",g,"' -A 1 ref_trna.fa > ref_gene.fasta",sep=""))
      g=gsub("\\(|\\)","",g)
    }
    
    
    vect<-c()
    
    for (sam in SAMPLES){
      
      system(paste("blastn -query /DATA/sdb1/Pierre/MitoZ/",sp,"/",sam,"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta -subject ref_gene.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send sequence' > tmp",sep=""))
      
      if (file.size("tmp") > 0){
        
        tmp<-read.table("tmp")
        
        if (length(unique(tmp$V1)) > 1){
          
          tmp=tmp[order(tmp$V3,decreasing = T),]
          
          system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",sam,"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_gene_indiv.fa",sep=""))
          system(paste("grep -w '",strsplit(as.character(tmp$V1[1]),";")[[1]][1],"' -A 1 ref_gene_indiv.fa > ref_gene_indiv_proper.fasta",sep=""))
          system(paste("blastn -query ref_gene_indiv_proper.fasta -subject ref_gene.fasta -outfmt '2 qseqid sseqid pident length qstart qend sstart send sequence' > tmp",sep=""))
          
          
        } else {
          
          system(paste("blastn -query /DATA/sdb1/Pierre/MitoZ/",sp,"/",sam,"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta -subject ref_gene.fasta -outfmt '2 qseqid sseqid pident length qstart qend sstart send sequence' > tmp",sep=""))
          
        }
        
        #system(paste("blastn -query /home/labosea1/MitoZ/ANALYSIS/",sp,"/Align/",sp,"_mitochondrial_align.fa -subject /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds -outfmt '2 qseqid sseqid pident length qstart qend sstart send sequence' > tmp",sep=""))
        tmp<-read.delim("tmp")
        
        data_weird=data.frame(QUERY=c(NA),
                              SUBJECT=c(NA),
                              SEQUENCE=c(NA))
        
        for (row in 1:nrow(tmp)){
          
          if (grepl("Query_",stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[1]) == T & grepl("Subject_",stri_remove_empty(strsplit(as.character(tmp[row+1,])," ")[[1]])[1]) == T){
            
            tmp1=seq(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[2],stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[4])
            tmp2=seq(stri_remove_empty(strsplit(as.character(tmp[row+1,])," ")[[1]])[2],stri_remove_empty(strsplit(as.character(tmp[row+1,])," ")[[1]])[4])  
            
            if (length(tmp2) > 1){
              if (((tmp2[2] - tmp2[1]) > 0) & ((tmp1[2] - tmp1[1]) < 0)){
                
                tmp1=rev(tmp1)
                tmp3=strsplit(stri_reverse(as.character(reverseComplement(DNAString(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[3])))),"")[[1]]
                reverse=1
                #tmp3=strsplit(stri_reverse(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[3]),"")[[1]]
                
              } else if (((tmp2[2] - tmp2[1]) < 0) & ((tmp1[2] - tmp1[1]) > 0)){
                
                tmp1=rev(tmp1)
                tmp3=strsplit(stri_reverse(as.character(reverseComplement(DNAString(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[3])))),"")[[1]]
                reverse=1
                #tmp3=strsplit(stri_reverse(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[3]),"")[[1]]
                
              } else {
                
                tmp3=strsplit(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[3],"")[[1]]
                reverse=0
              }
            } else {
              
              if (reverse==0){
                
                tmp3=strsplit(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[3],"")[[1]]
                
              } else {
                
                tmp3=as.character(reverseComplement(DNAString(stri_remove_empty(strsplit(as.character(tmp[row,])," ")[[1]])[3])))
                
              }
              
            }
            
            if (length(tmp1) != length(tmp2)){
              
              tmp1=tmp1[1:min(length(tmp1),length(tmp2))]
              tmp2=tmp2[1:min(length(tmp1),length(tmp2))]
              tmp3=tmp3[1:min(length(tmp1),length(tmp2))]
              
              
            }
            
            
            data_weird=rbind(data_weird,
                             data.frame(QUERY=tmp1,
                                        SUBJECT=tmp2,
                                        SEQUENCE=tmp3))
            
          }
          
        }
        
        data_weird=data_weird[-1,]
        data_weird=data_weird[order(data_weird$SUBJECT),]
        
        vect<-c(vect,
                paste(">",sam,sep=""),
                paste(
                  data_weird$SEQUENCE,
                  sep="",collapse="")
        )
        
      }
      
    }
    
    fileConn<-file(paste("GENE/",sp,"_",g,".fasta",sep=""))
    writeLines(vect,fileConn)
    close(fileConn)  
    
    system(paste("/usr/local/anaconda/bin/mafft GENE/",sp,"_",g,".fasta > GENE/",sp,"_",g,"_align.fa",sep=""))
    
    beeData<-read.multiFASTA(paste("GENE/",sp,"_",g,"_align.fa",sep=""))
    pdf(paste("GENE/",sp,"_",g,"align_data.pdf",sep=""),width=30,height=6)
    plot(beeData)
    dev.off()
    
    system(paste("/usr/local/anaconda/bin/trimal -in GENE/",sp,"_",g,"_align.fa -out GENE/",sp,"_",g,"_align_trim.fa -htmlout GENE/",sp,"_",g,"_align_trim.html -gt 1",sep=""))
    #system(paste("/usr/local/anaconda/bin/iqtree -s ",sp,"_mitochondrial_align_trim.fa -m MF -T 4",sep=""))
    system(paste("/usr/local/anaconda/bin/iqtree -s GENE/",sp,"_",g,"_align_trim.fa -m TN+F+I -T 1 -redo",sep=""))
    system(paste("/usr/local/anaconda/bin/raxmlHPC -m GTRGAMMA -s GENE/",sp,"_",g,"_align_trim.fa -n ",sp," -p 12345",sep=""))
    system(paste("mv RAxML_* GENE/.",sep=""))
    
    tree <- read.tree(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/GENE/",sp,"_",g,"_align_trim.fa.treefile",sep=""))
    image_info=data.frame(ind=tree$tip.label,
                          loc=substr(tree$tip.label,6,7),
                          sp=substr(tree$tip.label,1,5))
    
    image_info$loc=factor(image_info$loc,levels=c("Li","Mu","Fa","Ga"),
                          labels=c("Gulf of Lion",
                                   "Costa Calida",
                                   "Algarve",
                                   "Bay of Biscay")
    )
    p1<-ggtree(tree,col="black",size=1,alpha=0.5) %<+% image_info + 
      geom_tiplab(hjust=-0.15,size=3,fontface="bold",aes(color=loc)) +
      geom_tippoint(aes(color=loc),size=4,alpha=0.75) +
      theme_tree2() +
      scale_color_manual(name="Location",
                         values=RColorBrewer::brewer.pal(n=4,name="RdBu"))+
      scale_fill_manual(name="Species",values=wesanderson::wes_palette("Zissou1",n=length(levels(factor(image_info$sp))),type="continuous")) +
      #xlim(c(0,max(tree$edge.length))) +
      theme(legend.position = "bottom",
            plot.title=element_text(hjust=0.5)) +
      labs(title="iqtree")
    
    tree <- read.tree(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/GENE/RAxML_bestTree.",sp,sep=""))
    image_info=data.frame(ind=tree$tip.label,
                          loc=substr(tree$tip.label,6,7),
                          sp=substr(tree$tip.label,1,5))
    
    image_info$loc=factor(image_info$loc,levels=c("Li","Mu","Fa","Ga"),
                          labels=c("Gulf of Lion",
                                   "Costa Calida",
                                   "Algarve",
                                   "Bay of Biscay")
    )
    p2<-ggtree(tree,col="black",size=1,alpha=0.5) %<+% image_info + 
      geom_tiplab(hjust=-0.15,size=3,fontface="bold",aes(color=loc)) +
      geom_tippoint(aes(color=loc),size=4,alpha=0.75) +
      theme_tree2() +
      scale_color_manual(name="Location",
                         values=RColorBrewer::brewer.pal(n=4,name="RdBu"))+
      scale_fill_manual(name="Species",values=wesanderson::wes_palette("Zissou1",n=length(levels(factor(image_info$sp))),type="continuous")) +
      #xlim(c(0,max(tree$edge.length)+max(tree$edge.length)*0.15)) +
      theme(legend.position = "bottom",
            plot.title=element_text(hjust=0.5)) +
      labs(title="RAxML")
    
    pdf(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/GENE/",sp,"_",g,"_tree.pdf",sep=""),width=30,height=7.5)
    print(ggpubr::ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend = "bottom"))
    dev.off()
    
    pdf(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/GENE/",sp,"_",g,"_iqtree.pdf",sep=""),width=15,height=10)
    print(p1)
    dev.off()
    
    pdf(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/GENE/",sp,"_",g,"_tree_raxml.pdf",sep=""),width=15,height=10)
    print(p2)
    dev.off()
    
    system(paste("/usr/local/anaconda/bin/python3 /home/labosea1/MitoZ/ANALYSIS/run_mitoz.py ",sp," /home/labosea1/MitoZ/ANALYSIS/",sp,"/GENE/",sp,"_",g,"_align.fa ",g,sep=""))
    system(paste("rm VCF/vcf_",sp,".vcf",sep=""))
    
  }
  
  # dn/ds
  
  
  
  
}



#### ----

run_run =0

if (run_run == 1){
  for (sp in SPECIES){
    
    SAMPLES=list.files(paste("/DATA/sdb1/Pierre/MitoZ/",sp,sep=""),full.names=F,recursive=F)
    
    # Extract reference
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref.fa",sep=""))
    
    if (sp=="Cgale"){
      
      Ck141_593551
      
    } else if (sp=="Mmerl"){
      
      Ck141_672172
      
    } else if (sp=="Scant") {
      
      
    }
    
    setwd(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,sep=""))
    
    gene_file<-read.delim("/DATA/sdb1/Pierre/MitoZ/Cjuli/CjuliLi1/mitoz_annotate_PE.result/summary.txt")
    seq=strsplit(as.character(gene_file[1,])," ")[[1]][1]
    
    list_gene=c()
    type_gene=c()
    for (row in 1:nrow(gene_file)){
      if (seq == strsplit(as.character(gene_file[row,])," ")[[1]][1]){
        list_gene=c(list_gene,
                    stri_remove_empty(strsplit(as.character(gene_file[row,])," ")[[1]])[7])
        type_gene=c(type_gene,
                    stri_remove_empty(strsplit(as.character(gene_file[row,])," ")[[1]])[6])
      }
    }
    list_gene=list_gene[-1]
    type_gene=type_gene[-1]
    
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.cds /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_cds.fa",sep=""))
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.rrna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_rrna.fa",sep=""))
    system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",data_reference$REF[which(as.character(data_reference$SPECIES)==as.character(sp))],"/mitoz_annotate_PE.result/mitoz_annotate_PE.trna /home/labosea1/MitoZ/ANALYSIS/",sp,"/ref_trna.fa",sep=""))
    
    for (gene in 1:nrow(list_gene)){
      
      if (type_gene[gene]=="rRNA"){
        system(paste("grep -w '",list_gene[gene],"' -A 1 ref_rrna.fa > ref_gene.fasta",sep=""))
      } else if (type_gene[gene]=="tRNA"){
        system(paste("grep -w '",list_gene[gene],"' -A 1 ref_trna.fa > ref_gene.fasta",sep=""))
      } else {
        system(paste("grep -w '",list_gene[gene],"' -A 1 ref_cds.fa > ref_gene.fasta",sep=""))
      }
      
      vect<-c()
      for (sam in SAMPLES){
        
        system(paste("echo 'labosea1!' | sudo -S cp /DATA/sdb1/Pierre/MitoZ/",sp,"/",sam,"/mitoz_annotate_PE.result/mitoz_annotate_PE.fasta tmp_subject.fasta",sep=""))
        system(paste("blastn -query tmp_subject.fasta -subject ref_gene.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send' > tmp",sep=""))
        tmp<-read.table("tmp")
        file<-read.delim("tmp_subject.fasta",header=F)
        tmp=tmp[order(tmp$V4,decreasing = T),]
        
        if (tmp$V5[1]<tmp$V6[1]){
          sequence=substr(as.character(file[2,]),tmp$V5,tmp$V6)
        } else {
          print(sam)
          sequence=substr(as.character(file[2,]),tmp$V5,tmp$V6)
        }
        vect<-c(vect,
                paste(">",sam,sep=""),
                sequence)
      }
      
      fileConn<-file(paste(sp,"_",list_gene[gene],"_mitochondrial.fasta",sep=""))
      writeLines(vect,fileConn)
      close(fileConn)
      system(paste("/usr/local/anaconda/bin/mafft ",sp,"_",list_gene[gene],"_mitochondrial.fasta > ",sp,"_",list_gene[gene],"_mitochondrial_align",sep=""))
      
      beeData<-read.multiFASTA(paste(sp,"_",list_gene[gene],"_mitochondrial_align",sep=""))
      plot(beeData)
      
    }
    
  }
  
  
  
  
}
