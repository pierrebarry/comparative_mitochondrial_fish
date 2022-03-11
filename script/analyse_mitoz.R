## Compare species mitoz

library(apex)
library(Biostrings)
library(stringi)
library(ggtree)
library(tidyverse)
library(knitr)
library(kableExtra)
library(alignfigR)
library(ggrepel)
library(png)
library(grid)
library(apex)
library(Biostrings)
library(stringi)
library(ggtree)
library(tidyverse)
## Merge data frame

data_mitoch=data.frame(SP=c(NA),
                       GENE=c(NA),
                       LENGTH=c(NA),
                       BIALLELIC=c(NA),
                       POLYALLELIC=c(NA),
                       GAPS_ALL_INDIVIDUALS=c(NA),
                       GAPS=c(NA),
                       GC_CONTENT=c(NA),
                       PI_ALL=c(NA),
                       PI_LI=c(NA),
                       PI_MU=c(NA),
                       PI_FA=c(NA),
                       PI_GA=c(NA),
                       DXY_LI_MU=c(NA),
                       DXY_LI_FA=c(NA),
                       DXY_LI_GA=c(NA),
                       DXY_MU_FA=c(NA),
                       DXY_MU_GA=c(NA),
                       DXY_FA_GA=c(NA),		
                       DA_LI_MU=c(NA),
                       DA_LI_FA=c(NA),
                       DA_LI_GA=c(NA),
                       DA_MU_FA=c(NA),
                       DA_MU_GA=c(NA),
                       DA_FA_GA=c(NA),
                       FST_HUDSON_LI_MU=c(NA),
                       FST_HUDSON_LI_FA=c(NA),
                       FST_HUDSON_LI_GA=c(NA),
                       FST_HUDSON_MU_FA=c(NA),
                       FST_HUDSON_MU_GA=c(NA),
                       FST_HUDSON_FA_GA=c(NA),
                       FST_WEIRCOCKERHAM_LI_MU=c(NA),
                       FST_WEIRCOCKERHAM_LI_FA=c(NA),
                       FST_WEIRCOCKERHAM_LI_GA=c(NA),
                       FST_WEIRCOCKERHAM_MU_FA=c(NA),
                       FST_WEIRCOCKERHAM_MU_GA=c(NA),
                       FST_WEIRCOCKERHAM_FA_GA=c(NA),
                       TAJIMAD_ALL=c(NA),
                       TAJIMAD_LI=c(NA),
                       TAJIMAD_MU=c(NA),
                       TAJIMAD_FA=c(NA),
                       TAJIMAD_GA=c(NA)
)


SPECIES=c(
  "Aboye",
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
  "Styph"
)

for (sp in SPECIES){
  
  for (j in list.files(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/Pop_Gen",sep=""),full.names  = T)){
    
    whole<-read.table(j,header=T,sep=",")
    
    data_mitoch=rbind(data_mitoch,
                      whole)
    
  }
  
  
  
}

data_mitoch=data_mitoch[-1,]

type_gene=c()
for (i in 1:nrow(data_mitoch)){
  
  if (grepl("-rRNA",data_mitoch$GENE[i])){
    type_gene[i]="RNA"
  } else if (grepl("trn",data_mitoch$GENE[i])){
    type_gene[i]="RNA"
  } else {
    type_gene[i]="CDS"
  }
  
}

data_mitoch$type_gene=type_gene

for (i in 1:nrow(data_mitoch)){
  
  for (j in 26:37){
    
    if (is.na(data_mitoch[i,j])==T){
      data_mitoch[i,j]=0
    }
    
  }
  
}

data_mitoch$LENGTH=as.numeric(data_mitoch$LENGTH)
tt=data_mitoch[data_mitoch$GENE=="whole_mitochondrial",]
summary(aov(GAPS~SP, data=tt))

dd=c()

for (i in 1:nrow(data_mitoch[data_mitoch$GENE=="whole_mitochondrial",])){
  
  dd=c(dd,
       mean(as.numeric(data_mitoch[data_mitoch$GENE=="whole_mitochondrial",][i,10:13]),na.rm=T))
  
}

lfh<-readxl::read_excel("/home/labosea1/ARTICLE/mitoch/data/lfh.xlsx")
colnames(lfh)[3]="SP"
data_mitoch<-merge(data_mitoch,lfh[,c(3,5,6,7,8,9,10,11,12,13,14)],on=c("SP"))

pin_pis<-as.data.frame(readxl::read_excel("/home/labosea1/ARTICLE/mitoch/data/pin_pis.xlsx"))
ggplot(pin_pis,aes(x=SP,y=pin_pis,group=GENE)) +
  geom_boxplot() +
  geom_jitter(alpha=0.4,aes(col=GENE)) +
  geom_line() +
  theme_classic()
ggplot(pin_pis,aes(x=GENE,y=ratio_sites)) +
  geom_violin() +
  geom_jitter(alpha=0.4,aes(col=SP)) +
  #geom_line() +
  theme_classic() 
  #scale_y_log10()

data_mitoch=merge(data_mitoch,pin_pis,on=c("SP","GENE"))

ggplot(data_mitoch,aes(x=pin_pis,y=DXY_LI_GA)) +
  facet_wrap(SP~.,scales="free")+
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic()

ggplot(data_mitoch[data_mitoch$GENE!="whole_mitochondrial",],aes(x=pin_pis,y=PI_ALL)) +
  geom_point(aes(col=Parental_Care)) 





##

##

ggplot(data_mitoch[data_mitoch$type_gene=="CDS" & data_mitoch$GENE!="whole_mitochondrial",],aes(x=FST_HUDSON_LI_GA)) +
  geom_histogram() +
  facet_wrap(SP~.,scales="free")


tt<-data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",]
tt$pin_pis=as.numeric(tt$pin_pis)
tt$DXY_LI_GA=as.numeric(tt$DXY_LI_GA)
a<-lm(PI_ALL ~ ratio_sites, data=tt)
AIC(a)
a<-nlme::lme(FST_HUDSON_LI_GA ~ LENGTH, random = ~ 1 | SP/GAPS, data=tt[is.na(tt$pin_pis)==F,])
AIC(a)
summary(a)

a<-nlme::lme(FST_HUDSON_LI_GA ~ ratio_sites + pin_pis  + LENGTH, random = ~ 1 | SP/GAPS, data=tt[is.na(tt$pin_pis)==F & (tt$SP=="Gnige" | tt$SP=="Scant" | tt$SP=="Cgale" | tt$SP=="Cjuli" | tt$SP=="Dlabr" | tt$SP=="Lmorm"),])
summary(a)
a<-nlme::lme(PI_ALL ~ ratio_sites + pin_pis  + LENGTH, random = ~ 1 | SP/GAPS, data=tt[is.na(tt$pin_pis)==F,])
summary(a)

a<-nlme::lme(PI_ALL ~ ratio_sites + pin_pis, random = ~ 1 | SP/GAPS, data=tt[is.na(tt$pin_pis)==F,])
a<-nlme::lme(PI_ALL ~   pin_pis  + LENGTH, random = ~ 1 | SP/GAPS, data=tt[is.na(tt$pin_pis)==F,])

AIC(a)
summary(a)


save(data_mitoch,file="/home/labosea1/ARTICLE/mitoch/data/data_mitoch.Rdata")

##########


data_tmp=data.frame(SP=unique(data_mitoch$SP),
                    italic_species=c("italic('A. boyeri')",
                                     "italic('A. fallax')",
                                     "italic('C. galerita')",
                                     "italic('C. julis')",
                                     "italic('D. labrax')",
                                     "italic('D. puntazzo')",
                                     "italic('E. encrasicolus')",
                                     "italic('G. niger')",
                                     "italic('H. guttulatus')",
                                     "italic('L. budegassa')",
                                     "italic('L. mormyrus')",
                                     "italic('M. merluccius')",
                                     "italic('M. surmuletus')",
                                     "italic('P. erythrinus')",
                                     "italic('S. cabrilla')",
                                     "italic('S. cantharus')",
                                     "italic('S. cinereus')",
                                     "italic('S. pilchardus')",
                                     "italic('S. sarda')",
                                     "italic('S. typhle')")
)

data_mitoch=merge(data_mitoch,data_tmp,on=c("SP"))

tt<-data_mitoch[data_mitoch$GENE=="whole_mitochondrial",][,c(1,3,4,5,7,8)]
row.names(tt)=NULL
tt$SP=c("\\textit{A. boyeri}",
        "\\textit{A. fallax}",
        "\\textit{C. galerita}",
        "\\textit{C. julis}",
        "\\textit{D. labrax}",
        "\\textit{D. puntazzo}",
        "\\textit{E. encrasicolus}",
        "\\textit{G. niger}",
        "\textit{H. guttulatus}",
        "\\textit{L. budegassa}",
        "\\textit{L. mormyrus}",
        "\\textit{M. merluccius}",
        "\\textit{M. surmuletus}",
        "\\textit{P. erythrinus}",
        "\\textit{S. cabrilla}",
        "\\textit{S. cantharus}",
        "\\textit{S. cinereus}",
        "\\textit{S. pilchardus}",
        "\\textit{S. sarda}",
        "\\textit{S. typhle}")



print(knitr::kable(tt,"latex",escape=F),include.rownames=FALSE)

tt<-data_mitoch[data_mitoch$GENE=="whole_mitochondrial",][,c(1,9:13)]
for (i in seq(2,6)){
  tt[,i]=tt[,i]*100
}
tt$SP=c("\\textit{A. boyeri}",
        "\\textit{A. fallax}",
        "\\textit{C. galerita}",
        "\\textit{C. julis}",
        "\\textit{D. labrax}",
        "\\textit{D. puntazzo}",
        "\\textit{E. encrasicolus}",
        "\\textit{G. niger}",
        "\\textit{H. guttulatus}",
        "\\textit{L. budegassa}",
        "\\textit{L. mormyrus}",
        "\\textit{M. merluccius}",
        "\\textit{M. surmuletus}",
        "\\textit{P. erythrinus}",
        "\\textit{S. cabrilla}",
        "\\textit{S. cantharus}",
        "\\textit{S. cinereus}",
        "\\textit{S. pilchardus}",
        "\\textit{S. sarda}",
        "\\textit{S. typhle}")
print(knitr::kable(tt,"latex",escape=F),include.rownames=FALSE)

tt<-data_mitoch[data_mitoch$GENE=="whole_mitochondrial",][,c(1,26:31)]
#for (i in seq(2,6)){
#  tt[,i]=tt[,i]*100
#}
tt$SP=c("\\textit{A. boyeri}",
        "\\textit{A. fallax}",
        "\\textit{C. galerita}",
        "\\textit{C. julis}",
        "\\textit{D. labrax}",
        "\\textit{D. puntazzo}",
        "\\textit{E. encrasicolus}",
        "\\textit{G. niger}",
        "\\textit{H. guttulatus}",
        "\\textit{L. budegassa}",
        "\\textit{L. mormyrus}",
        "\\textit{M. merluccius}",
        "\\textit{M. surmuletus}",
        "\\textit{P. erythrinus}",
        "\\textit{S. cabrilla}",
        "\\textit{S. cantharus}",
        "\\textit{S. cinereus}",
        "\\textit{S. pilchardus}",
        "\\textit{S. sarda}",
        "\\textit{S. typhle}")
print(knitr::kable(tt,"latex",escape=F),include.rownames=FALSE)

tt<-data_mitoch[data_mitoch$GENE=="whole_mitochondrial",][,c(1,14:19)]
for (i in seq(2,6)){
  tt[,i]=tt[,i]*100
}
tt$SP=c("\\textit{A. boyeri}",
        "\\textit{A. fallax}",
        "\\textit{C. galerita}",
        "\\textit{C. julis}",
        "\\textit{D. labrax}",
        "\\textit{D. puntazzo}",
        "\\textit{E. encrasicolus}",
        "\\textit{G. niger}",
        "\\textit{H. guttulatus}",
        "\\textit{L. budegassa}",
        "\\textit{L. mormyrus}",
        "\\textit{M. merluccius}",
        "\\textit{M. surmuletus}",
        "\\textit{P. erythrinus}",
        "\\textit{S. cabrilla}",
        "\\textit{S. cantharus}",
        "\\textit{S. cinereus}",
        "\\textit{S. pilchardus}",
        "\\textit{S. sarda}",
        "\\textit{S. typhle}")
print(knitr::kable(tt,"latex",escape=F),include.rownames=FALSE)

tt<-data_mitoch[data_mitoch$GENE=="whole_mitochondrial",][,c(1,20:25)]
for (i in seq(2,6)){
  tt[,i]=tt[,i]*100
}
tt$SP=c("\\textit{A. boyeri}",
        "\\textit{A. fallax}",
        "\\textit{C. galerita}",
        "\\textit{C. julis}",
        "\\textit{D. labrax}",
        "\\textit{D. puntazzo}",
        "\\textit{E. encrasicolus}",
        "\\textit{G. niger}",
        "\\textit{H. guttulatus}",
        "\\textit{L. budegassa}",
        "\\textit{L. mormyrus}",
        "\\textit{M. merluccius}",
        "\\textit{M. surmuletus}",
        "\\textit{P. erythrinus}",
        "\\textit{S. cabrilla}",
        "\\textit{S. cantharus}",
        "\\textit{S. cinereus}",
        "\\textit{S. pilchardus}",
        "\\textit{S. sarda}",
        "\\textit{S. typhle}")
print(knitr::kable(tt,"latex",escape=F),include.rownames=FALSE)


tt<-data_mitoch[data_mitoch$GENE=="whole_mitochondrial",][,c(1,38:42)]
#for (i in seq(2,6)){
#  tt[,i]=tt[,i]*100
#}
tt$SP=c("\\textit{A. boyeri}",
        "\\textit{A. fallax}",
        "\\textit{C. galerita}",
        "\\textit{C. julis}",
        "\\textit{D. labrax}",
        "\\textit{D. puntazzo}",
        "\\textit{E. encrasicolus}",
        "\\textit{G. niger}",
        "\\textit{H. guttulatus}",
        "\\textit{L. budegassa}",
        "\\textit{L. mormyrus}",
        "\\textit{M. merluccius}",
        "\\textit{M. surmuletus}",
        "\\textit{P. erythrinus}",
        "\\textit{S. cabrilla}",
        "\\textit{S. cantharus}",
        "\\textit{S. cinereus}",
        "\\textit{S. pilchardus}",
        "\\textit{S. sarda}",
        "\\textit{S. typhle}")
print(knitr::kable(tt,"latex",escape=F),include.rownames=FALSE)


## length
data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(min(LENGTH),
            max(LENGTH),
            median(LENGTH),
            sd(LENGTH))

fastp<-read.table("/home/labosea1/ARTICLE/mitoch/data/Summary_fastp.txt",sep=";")
fastp=fastp[,c(1,16)]
colnames(fastp)=c("SAMPLES","GC")
fastp$SP=substr(fastp$SAMPLES,1,5)
tt<-tapply(fastp$GC,fastp$SP,median)

data_GC_nuclear=data.frame(SP=names(tt),
                           GC_nuclear=as.vector(tt))

data_GC=merge(data_GC_nuclear,
              data_mitoch[data_mitoch$GENE=="whole_mitochondrial",c("SP","GC_CONTENT","italic_species")],
              on=c("SP"))

colnames(data_GC)[3]="GC_mitoch"

pdf("/home/labosea1/ARTICLE/mitoch/figures/gc_nuclear_mitoch.pdf",width=7.5,height=7.5)
plot(
  ggplot(data_GC,aes(x=GC_mitoch,y=GC_nuclear)) +
    geom_point(alpha=0.4,size=3) +
    theme_classic() +
    xlab("Mitochondrial GC content (%)") +
    ylab("Nuclear GC content (%)") +
    geom_abline(slope=1,intercept=0,col='red',lty=2,alpha=0.5,lwd=1.5) +
    geom_text_repel(aes(label=italic_species),col="black",parse=T)
)
dev.off()

summary(lm(GC_nuclear~GC_mitoch,data=data_GC))


data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(min(GAPS),
            max(GAPS),
            median(GAPS),
            sd(GAPS))

## ALIGNMENT

#par(mfrow=c(5,4))
myplot=vector('list',20)
italic_sp=c("Atherina_boyeri",
            "Alosa fallax",
            "Coryphoblennius galerita",
            "Coris julis","
         Dicentrarchus labrax",
            "Diplodus puntazzo",
            "Engraulis encrasicolus",
            "Gobius niger",
            "Hippocampus guttulatus",
            "Lophius budegassa",
            "Lithognathus mormyrus",
            "Merluccius merluccius",
            "Mullus surmuletus",
            "Pagellus erythrinus",
            "Serranus cabrilla",
            "Spondyliosoma cantharus",
            "Symphodus cinereus",
            "Sardina pilchardus",
            "Sarda sarda",
            "Syngnathus typhle")
par(mfrow=c(4,5))
a=0
for (i in list.dirs("/home/labosea1/MitoZ/ANALYSIS/",recursive=F,full.names = F)){
  a=a+1
  beeData<-read.fasta(paste("/home/labosea1/MitoZ/ANALYSIS/",i,"/Align/",i,"_mitochondrial_align.fa",sep=""))
  colors= c("a" = "blue","c"="orange","g"="seagreen","t"="yellow","-"="black")
  #myplot[[i]]<-local({
  #plot.fasta(beeData,main=paste("expression(italic(",italic_sp[a],"))",sep=""))
  p1<-as.ggplot(plot(plot.fasta(beeData,main=expression(italic("A. boyeri")))))
  p2<-as.ggplot(plot(plot.fasta(beeData,main=expression(italic("C. julis")))))
  #})
}



data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  ggplot(aes(x=SP,y=LENGTH)) +
  geom_bar(stat='identity',fill="dodgerblue2",alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Mitochondrial scaffold length") +
  xlab("Species")

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>% 
  mutate(BIALLELIC = BIALLELIC / LENGTH) %>%
  ggplot(aes(x=SP,y=BIALLELIC)) +
  geom_bar(stat='identity',fill="dodgerblue2",alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Biallelic sites (%)") +
  xlab("Species")

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>% 
  mutate(POLYALLELIC = POLYALLELIC / LENGTH) %>%
  ggplot(aes(x=SP,y=POLYALLELIC)) +
  geom_bar(stat='identity',fill="dodgerblue2",alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Polyallelic sites (%)") +
  xlab("Species")

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>% 
  mutate(GAPS_ALL_INDIVIDUALS = GAPS_ALL_INDIVIDUALS / LENGTH) %>%
  ggplot(aes(x=SP,y=GAPS_ALL_INDIVIDUALS)) +
  geom_bar(stat='identity',fill="dodgerblue2",alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Gaps all individuals (%)") +
  xlab("Species")

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>% 
  mutate(GAPS = GAPS / LENGTH) %>%
  ggplot(aes(x=SP,y=GAPS)) +
  geom_bar(stat='identity',fill="dodgerblue2",alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Gaps (%)") +
  xlab("Species")

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>% 
  ggplot(aes(x=SP,y=GC_CONTENT)) +
  geom_bar(stat='identity',fill="dodgerblue2",alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("GC content (%)") +
  xlab("Species")

data_mitoch[,c(1,2,9,10,11,12,13)] %>%
  filter(GENE == "whole_mitochondrial") %>% 
  reshape2::melt(id.vars=c("SP","GENE")) %>%
  ggplot(aes(x=SP,y=value)) +
  geom_bar(position="dodge",stat='identity',aes(fill=variable),alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Genetic diversity (%)") +
  xlab("Species") +
  scale_fill_viridis_d(begin=0.1,end=0.9)

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  mutate(PI_ALL = PI_ALL * 100) %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(SP[which(min(PI_ALL)==PI_ALL)],
            min(PI_ALL),
            SP[which(max(PI_ALL)==PI_ALL)],
            max(PI_ALL),
            median(PI_ALL),
            sd(PI_ALL))

pdf("/home/labosea1/ARTICLE/mitoch/figures/pi_within_total.pdf",width=8.5,height=8.5)
plot(
  data_mitoch %>%
    filter(GENE == "whole_mitochondrial") %>%
    mutate(PI_MEAN = (PI_LI+PI_MU+PI_FA+PI_GA)*100/4) %>%
    mutate(PI_ALL = PI_ALL*100) %>%
    ggplot(aes(x=PI_ALL,y=PI_MEAN)) +
    geom_abline(slope=1,intercept=0,col='red',lty=2,alpha=0.5,lwd=1.5) +
    geom_point(alpha=0.5,size=4,aes(col=DXY_LI_GA)) +
    theme_classic() +
    xlab(expression(pi[total])) +
    ylab(expression(pi[within])) +
    scale_x_continuous(labels=seq(1,10),
                       breaks=seq(1,10)) +
    scale_y_continuous(labels=seq(1,10),
                       breaks=seq(1,10)) +
    scale_color_viridis_c(name=expression(d[XY]),begin=0.1,end=0.9) +
    geom_text_repel(aes(label=italic_species),col="black",parse=T,force=10)
  
)
dev.off()

data_mitoch[,c(1,2,14:19)] %>%
  filter(GENE == "whole_mitochondrial") %>% 
  reshape2::melt(id.vars=c("SP","GENE")) %>%
  ggplot(aes(x=SP,y=value)) +
  geom_bar(position="dodge",stat='identity',aes(fill=variable),alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Genetic divergence (%)") +
  xlab("Species") +
  scale_fill_viridis_d(begin=0.1,end=0.9)

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  mutate(DXY_LI_GA= DXY_LI_GA*100) %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(SP[which(min(DXY_LI_GA)==DXY_LI_GA)],
            min(DXY_LI_GA),
            SP[which(max(DXY_LI_GA)==DXY_LI_GA)],
            max(DXY_LI_GA),
            median(DXY_LI_GA),
            sd(DXY_LI_GA))

summary(lm(PI_ALL~DXY_LI_GA,data=data_mitoch))

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  filter(SP != "Aboye") %>%
  mutate(DXY_MU_FA= DXY_MU_FA*100) %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(SP[which(min(DXY_MU_FA,na.rm=T)==DXY_MU_FA)],
            min(DXY_MU_FA,na.rm=T),
            SP[which(max(DXY_MU_FA,na.rm=T)==DXY_MU_FA)],
            max(DXY_MU_FA,na.rm=T),
            median(DXY_MU_FA,na.rm=T),
            sd(DXY_MU_FA,na.rm=T))

summary(lm(DXY_MU_FA~DXY_LI_GA,data=data_mitoch))

plot(
  data_mitoch %>%
    filter(GENE == "whole_mitochondrial") %>%
    filter(SP !="Aboye") %>%
    mutate(DXY_LI_GA = DXY_LI_GA*100) %>%
    mutate(DXY_MU_FA= DXY_MU_FA*100) %>%
    ggplot(aes(y=DXY_LI_GA,x=DXY_MU_FA)) +
    geom_abline(slope=1,intercept=0,col='red',lty=2,alpha=0.5,lwd=1.5) +
    geom_point(alpha=0.5,size=4,col="chartreuse4") +
    theme_classic() +
    xlab(expression(pi[total])) +
    ylab(expression(pi[within])) +
    scale_x_continuous(labels=seq(1,10),
                       breaks=seq(1,10)) +
    scale_y_continuous(labels=seq(1,10),
                       breaks=seq(1,10)) +
    geom_text_repel(aes(label=italic_species),col="black",parse=T,force=10)
  
)

data_mitoch[,c(1,2,20:25)] %>%
  filter(GENE == "whole_mitochondrial") %>% 
  reshape2::melt(id.vars=c("SP","GENE")) %>%
  ggplot(aes(x=SP,y=value)) +
  geom_bar(position="dodge",stat='identity',aes(fill=variable),alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Absolute genetic divergence (%)") +
  xlab("Species") +
  scale_fill_viridis_d(begin=0.1,end=0.9)

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  mutate(DA_LI_GA= DA_LI_GA*100) %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(SP[which(min(DA_LI_GA)==DA_LI_GA)],
            min(DA_LI_GA),
            SP[which(max(DA_LI_GA)==DA_LI_GA)],
            max(DA_LI_GA),
            median(DA_LI_GA),
            sd(DA_LI_GA))

summary(lm(DXY_LI_GA~DA_LI_GA,data=data_mitoch))

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  filter(SP != "Aboye") %>%
  mutate(DA_MU_FA= DA_MU_FA*100) %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(SP[which(min(DA_MU_FA,na.rm=T)==DA_MU_FA)],
            min(DA_MU_FA,na.rm=T),
            SP[which(max(DA_MU_FA,na.rm=T)==DA_MU_FA)],
            max(DA_MU_FA,na.rm=T),
            median(DA_MU_FA,na.rm=T),
            sd(DA_MU_FA,na.rm=T))

summary(lm(DA_MU_FA~DA_LI_GA,data=data_mitoch))

data_mitoch$DXY_LI_GA=as.numeric(as.character(data_mitoch$DXY_LI_GA))

## PLot
data_mitoch_tmp=arrange(data_mitoch[data_mitoch$GENE=="whole_mitochondrial",],desc(DXY_LI_GA))
labels=as.character(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$italic_species)
p<-
  data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  filter(SP !="Aboye") %>%
  select("SP","DXY_LI_GA","DXY_MU_FA","italic_species") %>%
  mutate(SP=fct_reorder(SP,desc(DXY_LI_GA))) %>%
  mutate(DXY_LI_GA = DXY_LI_GA*100) %>%
  mutate(DXY_MU_FA= DXY_MU_FA*100) %>%
  reshape2::melt(id.vars=c("SP","italic_species")) %>%
  ggplot(aes(x=SP,y=value,fill=variable)) +
  geom_histogram(stat="identity",position=position_dodge(),alpha=0.5,col="grey") +
  theme_classic() +
  xlab("Species") +
  ylab(expression(d[XY])) +
  scale_y_continuous(expand=c(0,0),
                     labels=seq(1,10),
                     breaks=seq(1,10),
                     limits=c(0,8)) +
  scale_fill_manual(name=expression(d[XY]),
                    values=viridis::cividis(2,begin = 0.5,end=0.9),
                    labels=c("Between out","Between in")) +
  #scale_fill_viridis_d(begin=0.25,end=0.75) +
  scale_x_discrete(labels=parse(text=c(labels))) 

for (i in 1:nrow(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",])){
  img<-readPNG(paste("/home/labosea1/image/",data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$SP[i],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+
    annotation_custom(g,
                      xmin=i-0.5,
                      xmax=i+0.5,
                      ymin=max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DXY_LI_GA,na.rm=T)*100+max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DXY_LI_GA,na.rm=T)*100*0.01,
                      ymax=max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DXY_LI_GA,na.rm=T)*100+max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DXY_LI_GA,na.rm=T)*100*0.1)
}

p_dxy<-p
pdf("/home/labosea1/ARTICLE/mitoch/figures/dxy_mitoch_continuum.pdf",width=17.5,height=5)
print(p_dxy)
dev.off()

data_mitoch_tmp=arrange(data_mitoch[data_mitoch$GENE=="whole_mitochondrial",],desc(DA_LI_GA))
labels=as.character(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$italic_species)
p<-
  data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  filter(SP !="Aboye") %>%
  select("SP","DA_LI_GA","DA_MU_FA","italic_species") %>%
  mutate(SP=fct_reorder(SP,desc(DA_LI_GA))) %>%
  mutate(DA_LI_GA = DA_LI_GA*100) %>%
  mutate(DA_MU_FA= DA_MU_FA*100) %>%
  reshape2::melt(id.vars=c("SP","italic_species")) %>%
  ggplot(aes(x=SP,y=value,fill=variable)) +
  geom_histogram(stat="identity",position=position_dodge(),alpha=0.5,col="grey") +
  theme_classic() +
  xlab("Species") +
  ylab(expression(d[a])) +
  scale_y_continuous(expand=c(0,0),
                     labels=seq(1,10),
                     breaks=seq(1,10),
                     limits=c(0,8)) +
  scale_fill_manual(name=expression(d[a]),
                    values=viridis::cividis(2,begin = 0.5,end=0.9),
                    labels=c("Between out","Between in")) +
  #scale_fill_viridis_d(begin=0.25,end=0.75) +
  scale_x_discrete(labels=parse(text=c(labels))) 

for (i in 1:nrow(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",])){
  img<-readPNG(paste("/home/labosea1/image/",data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$SP[i],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+
    annotation_custom(g,
                      xmin=i-0.5,
                      xmax=i+0.5,
                      ymin=max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DA_LI_GA,na.rm=T)*100+max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DA_LI_GA,na.rm=T)*100*0.01,
                      ymax=max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DA_LI_GA,na.rm=T)*100+max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$DA_LI_GA,na.rm=T)*100*0.1)
}

p_da<-p
pdf("/home/labosea1/ARTICLE/mitoch/figures/da_mitoch_continuum.pdf",width=17.5,height=5)
print(p_da)
dev.off()

data_mitoch_tmp=arrange(data_mitoch[data_mitoch$GENE=="whole_mitochondrial",],desc(FST_HUDSON_LI_GA))
labels=as.character(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$italic_species)
p<-
  data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  filter(SP !="Aboye") %>%
  select("SP","FST_HUDSON_LI_GA","FST_HUDSON_MU_FA","italic_species") %>%
  mutate(SP=fct_reorder(SP,desc(FST_HUDSON_LI_GA))) %>%
  mutate(FST_HUDSON_LI_GA = ifelse(FST_HUDSON_LI_GA<0,0,FST_HUDSON_LI_GA)) %>%
  reshape2::melt(id.vars=c("SP","italic_species")) %>%
  ggplot(aes(x=SP,y=value,fill=variable)) +
  geom_histogram(stat="identity",position=position_dodge(),alpha=0.5,col="grey") +
  theme_classic() +
  xlab("Species") +
  ylab(expression(F[ST])) +
  scale_y_continuous(expand=c(0,0),
                     labels=seq(0,1,by=0.1),
                     breaks=seq(0,1,by=0.1),
                     limits=c(0,1.15)) +
  scale_fill_manual(name=expression(F[ST]),
                    values=viridis::cividis(2,begin = 0.5,end=0.9),
                    labels=c("Between out","Between in")) +
  #scale_fill_viridis_d(begin=0.25,end=0.75) +
  scale_x_discrete(labels=parse(text=c(labels))) 

for (i in 1:nrow(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",])){
  img<-readPNG(paste("/home/labosea1/image/",data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$SP[i],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p<-p+
    annotation_custom(g,
                      xmin=i-0.5,
                      xmax=i+0.5,
                      ymin=max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$FST_HUDSON_LI_GA,na.rm=T)+max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$FST_HUDSON_LI_GA,na.rm=T)*0.01,
                      ymax=max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$FST_HUDSON_LI_GA,na.rm=T)+max(data_mitoch_tmp[data_mitoch_tmp$SP!="Aboye",]$FST_HUDSON_LI_GA,na.rm=T)*0.1)
}

p_fst<-p
pdf("/home/labosea1/ARTICLE/mitoch/figures/fst_mitoch_continuum.pdf",width=17.5,height=5)
print(p_fst)
dev.off()

pdf("/home/labosea1/ARTICLE/mitoch/figures/fst_dxy_da_mitoch_continuum.pdf",width=17.5,height=15)
print(ggpubr::ggarrange(p_fst,
                        p_dxy,
                        p_da,
                        nrow=3))
dev.off()

data_mitoch[,c(1,2,26:31)] %>%
  filter(GENE == "whole_mitochondrial") %>% 
  reshape2::melt(id.vars=c("SP","GENE")) %>%
  ggplot(aes(x=SP,y=value)) +
  geom_bar(position="dodge",stat='identity',aes(fill=variable),alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,1)) +
  ylab("Hudson's Fst") +
  xlab("Species") +
  scale_fill_viridis_d(begin=0.1,end=0.9) 


data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  filter(FST_HUDSON_LI_GA>0) %>%
  #filter(SP != "Aboye") %>%
  mutate(FST_HUDSON_LI_GA= FST_HUDSON_LI_GA) %>%
  #select("SP","LENGTH") %>%
  #kable() %>%
  summarize(
    min(FST_HUDSON_LI_GA[FST_HUDSON_LI_GA>0],na.rm=T),
    SP[which(max(FST_HUDSON_LI_GA,na.rm=T)==FST_HUDSON_LI_GA)],
    max(FST_HUDSON_LI_GA,na.rm=T),
    median(FST_HUDSON_LI_GA,na.rm=T),
    sd(FST_HUDSON_LI_GA,na.rm=T))

summary(lm(DXY_LI_GA~FST_HUDSON_LI_GA,data=data_mitoch))
summary(lm(DA_LI_GA~FST_HUDSON_LI_GA,data=data_mitoch))

data_mitoch[,c(1,2,38:42)] %>%
  filter(GENE == "whole_mitochondrial") %>% 
  reshape2::melt(id.vars=c("SP","GENE")) %>%
  ggplot(aes(x=SP,y=value)) +
  geom_bar(position="dodge",stat='identity',aes(fill=variable),alpha=0.4) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  ylab("Tajima's D") +
  xlab("Species") +
  scale_fill_viridis_d(begin=0.1,end=0.9) 

data_mitoch %>%
  filter(GENE == "whole_mitochondrial") %>%
  select("SP","TAJIMAD_ALL")

###

tmp<-data_mitoch[data_mitoch$GENE=="whole_mitochondrial",]
load(file="/home/labosea1/ARTICLE/mitoch/data/dxy_all.Rdata")
dxy_data<-merge(tmp,dxy_data,on=c("SP"))
load(file="/home/labosea1/ARTICLE/mitoch/data/Summary_GenomeScope.Rdata")
dxy_data$pi_nuc=tapply(as.numeric(Summary_GenomeScope$Heterozygosity),Summary_GenomeScope$Species,median)[-c(2,7,8,14)]
fst<-as.data.frame(readxl::read_excel("/home/labosea1/ARTICLE/mitoch/data/FST.xlsx"))
colnames(fst)[2]="SP"
dxy_data<-merge(dxy_data,fst,on=c("SP"))

dxy_data$italic_species=c("italic('A. boyeri')",
                                                             #"italic('A. fallax')",
                                                             "italic('C. galerita')",
                                                             "italic('C. julis')",
                                                             "italic('D. labrax')",
                                                             "italic('D. puntazzo')",
                                                             #"italic('E. encrasicolus')",
                                                             #"italic('G. niger')",
                                                             "italic('H. guttulatus')",
                                                             "italic('L. budegassa')",
                                                             "italic('L. mormyrus')",
                                                             "italic('M. merluccius')",
                                                             "italic('M. surmuletus')",
                                                             "italic('P. erythrinus')",
                                                             "italic('S. cabrilla')",
                                                             "italic('S. cantharus')",
                                                             "italic('S. cinereus')",
                                                             "italic('S. pilchardus')",
                                                             "italic('S. sarda')",
                                                             "italic('S. typhle')")

dxy_data$dxy_Li_Ga=as.numeric(dxy_data$dxy_Li_Ga)
dxy_data$da_Li_Ga=as.numeric(dxy_data$da_Li_Ga)

p_dxy<-plot(
  ggplot(dxy_data,aes(y=DXY_LI_GA*100,x=dxy_Li_Ga)) +
    geom_point(alpha=0.4,size=3) +
    theme_classic() +
    xlab("Nuclear dxy (%)") +
    ylab("Mitochondrial dxy (%)") +
    geom_abline(slope=1,intercept=0,col='red',lty=2,alpha=0.5,lwd=1.5) +
    geom_text_repel(aes(label=italic_species),col="black",parse=T) +
    labs(title=expression(d[XY])) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour="black",fill=NA))
)

p_da<-plot(
  ggplot(dxy_data,aes(y=DA_LI_GA*100,x=da_Li_Ga)) +
    geom_point(alpha=0.4,size=3) +
    theme_classic() +
    xlab("Nuclear da (%)") +
    ylab("Mitochondrial da (%)") +
    geom_abline(slope=1,intercept=0,col='red',lty=2,alpha=0.5,lwd=1.5) +
    geom_text_repel(aes(label=italic_species),col="black",parse=T)+
    labs(title=expression(d[a])) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour="black",fill=NA))
)

p_fst<-plot(
  ggplot(dxy_data,aes(y=FST_HUDSON_LI_GA,x=FST_LI_GA_WEIGHTED)) +
    geom_point(alpha=0.4,size=3) +
    theme_classic() +
    xlab("Nuclear Fst") +
    ylab("Mitochondrial Fst") +
    geom_abline(slope=1,intercept=0,col='red',lty=2,alpha=0.5,lwd=1.5) +
    geom_text_repel(aes(label=italic_species),col="black",parse=T)+
    labs(title=expression(F[ST])) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour="black",fill=NA))
)

p_pi<-plot(
  ggplot(dxy_data,aes(y=PI_ALL*100,x=pi_nuc)) +
    geom_point(alpha=0.4,size=3) +
    theme_classic() +
    xlab("Nuclear pi (%)") +
    ylab("Mitochondrial pi (%)") +
    geom_abline(slope=1,intercept=0,col='red',lty=2,alpha=0.5,lwd=1.5) +
    geom_text_repel(aes(label=italic_species),col="black",parse=T)+
    labs(title=expression(pi)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour="black",fill=NA))
)

pdf("/home/labosea1/ARTICLE/mitoch/figures/fst_dxy_da_pi.pdf",width=12.5,height=12.5)
ggpubr::ggarrange(p_pi,
                  p_fst,
                  p_dxy,
                  p_da,
                  nrow=2,
                  ncol=2,
                  labels = c("A","B","C","D"))
dev.off()

summary(lm(PI_ALL~pi_nuc,data=dxy_data))
summary(lm(FST_HUDSON_LI_GA~FST_LI_GA_WEIGHTED,data=dxy_data))
summary(lm(DA_LI_GA~da_Li_Ga,data=dxy_data))
summary(lm(DXY_LI_GA~dxy_Li_Ga,data=dxy_data))

##

coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Aboye",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Afall",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Cjuli",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Cgale",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Dlabr",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Dpunt",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Eencr",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Gnige",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Hgutt",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Lbude",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Lmorm",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Mmerl",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Msurm",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Peryt",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scabr",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scant",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scine",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Spilc",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Ssard",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Styph",])))[2,4]


coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Aboye",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Afall",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Cjuli",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Cgale",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Dlabr",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Dpunt",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Eencr",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Gnige",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Hgutt",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Lbude",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Lmorm",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Mmerl",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Msurm",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Peryt",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scabr",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scant",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scine",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Spilc",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Ssard",])))[2,4]
coefficients(summary(lm(DXY_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Styph",])))[2,4]

coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Aboye",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Afall",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Cjuli",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Cgale",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Dlabr",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Dpunt",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Eencr",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Gnige",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Hgutt",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Lbude",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Lmorm",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Mmerl",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Msurm",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Peryt",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scabr",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scant",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Scine",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Spilc",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Ssard",])))[2,4]
coefficients(summary(lm(FST_HUDSON_LI_GA~type_gene,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$SP=="Styph",])))[2,4]





summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Aboye",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Afall",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Cgale",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Cjuli",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Dlabr",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Dpunt",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Eencr",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Gnige",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Hgutt",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Lbude",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Lmorm",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Mmerl",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Msurm",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Peryt",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Scabr",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Scant",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Scine",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Spilc",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Ssard",]))
summary(lm(FST_HUDSON_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP=="Styph",]))

library(scmamp)
data(data_gh_2008)

for (i in 1:nrow(data_mitoch)){
 
  for (j in 28){
    
    data_mitoch[i,j]=as.numeric(data_mitoch[i,j])/as.numeric(data_mitoch[data_mitoch$GENE=="whole_mitochondrial" & data_mitoch$SP==data_mitoch$SP[i],"FST_HUDSON_LI_GA"])
    
  }
  
}

ggplot(data_mitoch[data_mitoch$type_gene=="CDS",],aes(x=GENE,y=DXY_LI_GA,group=SP)) +
  geom_point(alpha=0.4,aes(col=SP)) +
  #geom_line(aes(col=SP)) +
  geom_violin() +
  theme_classic() +
  ylim(c(0,2))

for (i in head(unique(levels(factor(data_mitoch$SP))),-1)){
  
  for (j in unique(levels(factor(data_mitoch$SP)))[-1]){
    
    data_tmp1=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP==i,]
    data_tmp1=data_tmp1[order(data_tmp1$DXY_LI_GA),]
    data_tmp1$order=seq(1,nrow(data_tmp1))
    data_tmp2=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS" & data_mitoch$SP==j,]
    data_tmp2=data_tmp2[order(data_tmp2$DXY_LI_GA),]
    data_tmp2$order=seq(1,nrow(data_tmp2))
    
    data_tmp=rbind(data_tmp1[,c("SP","GENE","order")],
                   data_tmp2[,c("SP","GENE","order")])
    
    data_tmp$SP=factor(data_tmp$SP)
    data_tmp$GENE=factor(data_tmp$GENE)
    data_tmp$order=factor(data_tmp$order)
    #m1<-nlme::lme(order ~ GENE, random = ~ GENE | SP, data=data_tmp)
    #m2<-nlme::lme(order ~ GENE, random = ~ 1 | SP, data=data_tmp)
    #anova(m1,m2)
    ggplot(data_tmp,aes(x=GENE,y=order,group=SP)) +
      geom_point(aes(col=SP)) +
      geom_line(aes(col=SP)) +
      theme_classic()
    
    summary(lm(order~GENE:SP,data=data_tmp))
    
    readline(paste(i,"-",j,sep=""))
    
  }
  
}


cor_data=matrix(rep(0,20*20),nrow=20,ncol=20)
row.names(cor_data)=levels(factor(data_mitoch$SP))
colnames(cor_data)=levels(factor(data_mitoch$SP))

for (i in unique(levels(factor(data_mitoch$SP)))){
  
  for (j in unique(levels(factor(data_mitoch$SP)))){
    
    data_tmp1=data_mitoch[data_mitoch$GENE!="whole_mitochondrial"  & data_mitoch$type_gene=="CDS" & data_mitoch$SP==i,]
    data_tmp2=data_mitoch[data_mitoch$GENE!="whole_mitochondrial"  & data_mitoch$type_gene=="CDS" & data_mitoch$SP==j,]

    data_tmp1=data_tmp1[,c("GENE","PI_ALL")]
    data_tmp2=data_tmp2[,c("GENE","PI_ALL")]
    colnames(data_tmp1)=c("GENE","dxy1")
    colnames(data_tmp2)=c("GENE","dxy2")
    dd=merge(data_tmp1,
             data_tmp2,
             on=c("GENE"))
    
    if (nrow(dd)>=1){
      #print(ggplot(dd,aes(x=dxy1,y=dxy2)) +
      #  geom_point() +
      #  theme_classic()
      #)
      a<-summary(lm(dxy2~dxy1,data=dd))
      cor_data[which(row.names(cor_data)==i),which(colnames(cor_data)==j)]=coefficients(a)[2,4]
      #readline(paste(i,"-",j,sep=""))
    }

    
  }
  
}

options(digits=3)

#for (i in 1:nrow(cor_data)){
#  cor_data[,i]=round(cor_data[,i],3)
#}

cor_data=format(cor_data,scientific=T)
for (i in 1:nrow(cor_data)){
  cor_data[,i]=paste("\\tiny{${",cor_data[,i],"}$}",sep="")
}

for (i in 1:nrow(cor_data)){
  for (j in 1:ncol(cor_data)){
    if (j>=i){
      cor_data[i,j]="-"
    }
  }
}

#xtable(cor_data)

knitr::kable(cor_data,"latex",escape=F)






dd=data.frame(SP=c(NA),
                PVALUE=c(NA))
for (i in 1:nrow(cor_data)){
  
  dd=rbind(dd,
           data.frame(SP=row.names(cor_data)[i],
                      PVALUE=cor_data[i,]))
  print(length(cor_data[i,][cor_data[i,]<=0.05]))
  
}

dd=dd[-1,]
dd$PVALUE=log(dd$PVALUE)
ggplot(dd,aes(x=PVALUE)) +
  geom_histogram() +
  theme_classic() +
  facet_wrap(SP ~ .) +
  geom_vline(xintercept = log(0.05))

# RAtio 
data_mitoch_ratio<-data_mitoch

for (i in 1:nrow(data_mitoch)){
  
  for (j in seq(9,37)){
    
    data_mitoch_ratio[i,j]=data_mitoch[i,j]/(data_mitoch[data_mitoch$GENE=="whole_mitochondrial" & data_mitoch$SP==data_mitoch$SP[i],j])
    
  }
  
}


#pi

data_mitoch[,c(1,2,3,7,9,43)] %>%
  filter(GENE != "whole_mitochondrial") %>% 
  filter((GAPS/LENGTH)<0.1) %>%
  #filter(SP=="Scant") %>%
  #filter(type_gene=="CDS") %>%
  reshape2::melt(id.vars=c("SP","GENE","LENGTH","GAPS","type_gene")) %>%
  ggplot(aes(x=GENE,y=value,group=SP)) +
  #geom_point(aes(col=SP),alpha=0.4) +
  geom_line(aes(col=SP)) +
  #facet_wrap(SP~.,nrow=20) +
  theme_classic() +
  ylab(expression(pi)) +
  xlab("Gene") +
  scale_color_viridis_d(begin=0.1,end=0.9) 

#dxy

data_mitoch[,c(1,2,3,7,16,43)] %>%
  filter(GENE != "whole_mitochondrial") %>% 
  filter((GAPS/LENGTH)<0.1) %>%
  #filter(SP=="Scant") %>%
  #filter(type_gene=="CDS") %>%
  reshape2::melt(id.vars=c("SP","GENE","LENGTH","GAPS","type_gene")) %>%
  ggplot(aes(x=GENE,y=value,group=SP)) +
  #geom_point(aes(col=SP),alpha=0.4) +
  geom_line(aes(col=SP)) +
  #facet_wrap(SP~.,nrow=20) +
  theme_classic() +
  ylab(expression(d[XY])) +
  xlab("Gene") +
  scale_color_viridis_d(begin=0.1,end=0.9) 

#da

data_mitoch[,c(1,2,3,7,22,43)] %>%
  filter(GENE != "whole_mitochondrial") %>% 
  filter((GAPS/LENGTH)<0.1) %>%
  #filter(SP=="Scant") %>%
  #filter(type_gene=="CDS") %>%
  reshape2::melt(id.vars=c("SP","GENE","LENGTH","GAPS","type_gene")) %>%
  ggplot(aes(x=GENE,y=value,group=SP)) +
  #geom_point(aes(col=SP),alpha=0.4) +
  geom_line(aes(col=SP)) +
  #facet_wrap(SP~.,nrow=20) +
  theme_classic() +
  ylab(expression(d[a])) +
  xlab("Gene") +
  scale_color_viridis_d(begin=0.1,end=0.9) 


#fst

data_mitoch[,c(1,2,3,7,28,43)] %>%
  filter(GENE != "whole_mitochondrial") %>% 
  filter((GAPS/LENGTH)<0.1) %>%
  #filter(SP=="Scant") %>%
  #filter(type_gene=="CDS") %>%
  reshape2::melt(id.vars=c("SP","GENE","LENGTH","GAPS","type_gene")) %>%
  ggplot(aes(x=GENE,y=value,group=SP)) +
  #geom_point(aes(col=SP),alpha=0.4) +
  geom_line(aes(col=SP)) +
  #facet_wrap(SP~.,nrow=20) +
  theme_classic() +
  ylab(expression(F[ST])) +
  xlab("Gene") +
  scale_color_viridis_d(begin=0.1,end=0.9) 




for (i in levels(factor(tt$SP))){
  tt_1<-tt[tt$SP!=i,]
  a<-nlme::lme(DXY_LI_GA ~ LENGTH , random = ~ 1 | SP/GAPS/GENE, data=tt_1)
  #summary(a)
  if (coefficients(summary(a))[2,5]<0.05){
    print(coefficients(summary(a))[2,5])
  }
}

tt<-data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",]

a<-nlme::lme(FST_HUDSON_LI_GA ~ LENGTH, random = ~ 1 | SP, data=tt[tt$SP!="Mmerl" & tt$SP!="Scant",])
AIC(a)
summary(a)

a<-nlme::lme(FST_HUDSON_LI_GA ~ LENGTH, random = ~ 1 | GAPS, data=tt[tt$SP!="Mmerl" & tt$SP!="Scant",])
AIC(a)
summary(a)

a<-nlme::lme(FST_HUDSON_LI_GA ~ LENGTH, random = ~ 1 | SP/GAPS, data=tt[tt$SP!="Mmerl" & tt$SP!="Scant",])
AIC(a)
summary(a)

a<-nlme::lme(FST_HUDSON_LI_GA ~ LENGTH, random = ~ 1 | SP, data=tt)
AIC(a)
summary(a)

for (i in levels(factor(tt$SP))){
  tt_1<-tt[tt$SP!=i,]
  a<-nlme::lme(DXY_LI_GA ~ LENGTH , random = ~ 1 | SP/GAPS/GENE, data=tt_1)
#summary(a)
  if (coefficients(summary(a))[2,5]<0.05){
    print(coefficients(summary(a))[2,5])
  }
}


all_comb<- combn(length(levels(factor(tt$SP))),2)

for (i in 1:ncol(all_comb)){
  tt_1<-tt[tt$SP!=levels(factor(tt$SP))[all_comb[1,i]] & tt$SP!=levels(factor(tt$SP))[all_comb[2,i]],]
  a<-nlme::lme(DXY_LI_GA ~ LENGTH , random = ~ 1 | SP/GAPS/GENE, data=tt_1)
  #summary(a)
  if (coefficients(summary(a))[2,5]<0.05){
    print(coefficients(summary(a))[2,5])
    
  }
}

all_comb<- combn(length(levels(factor(tt$SP))),3)

for (i in 1:ncol(all_comb)){
  tt_1<-tt[tt$SP!=levels(factor(tt$SP))[all_comb[1,i]] & tt$SP!=levels(factor(tt$SP))[all_comb[2,i]] &  tt$SP!=levels(factor(tt$SP))[all_comb[3,i]],]
  a<-nlme::lme(DXY_LI_GA ~ LENGTH , random = ~ 1 | SP/GAPS/GENE, data=tt_1)
  #summary(a)
  if (coefficients(summary(a))[2,5]<0.05){
    print(coefficients(summary(a))[2,5])
    print(all_comb[,i])
    
  }
}

a<-lm(DA_LI_GA~LENGTH,data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",])
AIC(a)
summary(a)
a<-nlme::lme(DXY_LI_GA ~ LENGTH , random = ~ 1 | SP/GAPS/GENE, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",])
AIC(a)
summary(a)
a<-nlme::lme(DA_LI_GA ~ LENGTH , random = ~ SP | SP, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",])
AIC(a)
summary(a)


dxy_li_ga_ref=c()
da_li_ga_ref=c()
pi_li_ga_ref=c()
fst_li_ga_ref=c()


for (i in 1:nrow(data_mitoch)){
  
  if (data_mitoch$GENE[i]=="whole_mitochondrial"){
    dxy_li_ga_ref=c(dxy_li_ga_ref,NA)
    da_li_ga_ref=c(da_li_ga_ref,NA)
    pi_li_ga_ref=c(pi_li_ga_ref,NA)
    fst_li_ga_ref=c(fst_li_ga_ref,NA)
  } else {
    dxy_li_ga_ref=c(dxy_li_ga_ref,data_mitoch[data_mitoch$GENE=="whole_mitochondrial" & data_mitoch$SP==data_mitoch$SP[i],]$DXY_LI_GA)
    da_li_ga_ref=c(da_li_ga_ref,data_mitoch[data_mitoch$GENE=="whole_mitochondrial" & data_mitoch$SP==data_mitoch$SP[i],]$DA_LI_GA)
    pi_li_ga_ref=c(pi_li_ga_ref,data_mitoch[data_mitoch$GENE=="whole_mitochondrial" & data_mitoch$SP==data_mitoch$SP[i],]$PI_ALL)
    fst_li_ga_ref=c(fst_li_ga_ref,data_mitoch[data_mitoch$GENE=="whole_mitochondrial" & data_mitoch$SP==data_mitoch$SP[i],]$FST_HUDSON_LI_GA)
  }
  
}

data_mitoch$dxy_li_ga_ref=dxy_li_ga_ref
data_mitoch$da_li_ga_ref=da_li_ga_ref
data_mitoch$pi_li_ga_ref=pi_li_ga_ref
data_mitoch$fst_li_ga_ref=fst_li_ga_ref

a<-nlme::lme(DA_LI_GA ~ LENGTH , random = ~ 1 | da_li_ga_ref, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",])
AIC(a)
summary(a)
a<-nlme::lme(DXY_LI_GA ~ LENGTH + dxy_li_ga_ref, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",])
AIC(a)
summary(a)

summary(lm(DXY_LI_GA ~ LENGTH + dxy_li_ga_ref, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",]))
summary(lm(DA_LI_GA ~ LENGTH + da_li_ga_ref, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",]))
summary(lm(PI_ALL ~ LENGTH + pi_li_ga_ref, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",]))
summary(lm(FST_HUDSON_LI_GA ~ LENGTH + fst_li_ga_ref, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial" & data_mitoch$type_gene=="CDS",]))

a<-nlme::lme(DA_LI_GA ~ -1+GENE, random = ~ 1 | SP/GAPS, data=data_mitoch[data_mitoch$GENE!="whole_mitochondrial",])
summary(a)

dd=data.frame(GENE=substr(row.names(summary(a)$tTable),5,nchar(row.names(summary(a)$tTable))),
              ESTIMATE=summary(a)$tTable[,2],
              SD=summary(a)$tTable[,1])

type_gene=c()
for (i in 1:nrow(dd)){
  if ("trn" %in% as.character(dd$GENE[i])){
    type_gene[i]="RNA"
  } else if("rRNA" %in% as.character(dd$GENE[i])){
    type_gene[i]="RNA"
  } else{
    type_gene[i]="CDS"
  }
}

dd$type_gene=type_gene

mutate(dd,GENE=fct_reorder(GENE,ESTIMATE)) %>%
  ggplot(aes(x=GENE,y=ESTIMATE)) +
    geom_pointrange(aes(ymin=ESTIMATE-SD^(2),ymax=ESTIMATE+SD^(2)),aes(col=type_gene)) +
    theme_classic()

SPECIES=c("Aboye",
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
          "Styph")
myplots<-vector('list',20)
a=0
for (sp in SPECIES){
  
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
  a=a+1
  myplots[[a]] <- local({
    p1<-ggtree(tree,col="black",size=1,alpha=0.5) %<+% image_info + 
      #geom_tiplab(hjust=-0.15,size=3,fontface="bold",aes(color=loc)) +
      geom_tippoint(aes(color=loc),size=6,alpha=0.75) +
      theme_tree2() +
      scale_color_manual(name="Location",
                         values=RColorBrewer::brewer.pal(n=4,name="RdBu"))+
      scale_fill_manual(name="Species",values=wesanderson::wes_palette("Zissou1",n=length(levels(factor(image_info$sp))),type="continuous")) +
      #xlim(c(0,max(tree$edge.length))) +
      theme(legend.position = "bottom",
            plot.title=element_text(hjust=0.5)) 
      #labs(title="iqtree")
    p1
    img<-png::readPNG(paste("/home/labosea1/image/",sp,".png",sep=""))
    g <- grid::rasterGrob(img, interpolate=TRUE)
    
    data_fake=data.frame(X=seq(0,10,by=1),
                         Y=seq(0,10,by=1))
    p<-ggplot(data_fake,aes(x=X,y=Y)) +
      theme_void() +
      annotation_custom(g) +
      ylim(c(25,100))
    #p
    print(ggpubr::ggarrange(p,p1,nrow=2,heights = c(0.1,0.95)))
  })
  p1<-ggtree(tree,col="black",size=1,alpha=0.5) %<+% image_info + 
    #geom_tiplab(hjust=-0.15,size=3,fontface="bold",aes(color=loc)) +
    geom_tippoint(aes(color=loc),size=6,alpha=0.75) +
    theme_tree2() +
    scale_color_manual(name="Location",
                       values=RColorBrewer::brewer.pal(n=4,name="RdBu"))+
    scale_fill_manual(name="Species",values=wesanderson::wes_palette("Zissou1",n=length(levels(factor(image_info$sp))),type="continuous")) +
    #xlim(c(0,max(tree$edge.length))) +
    theme(legend.position = "none",
          plot.title=element_text(hjust=0.5)) 
  #labs(title="iqtree")
  p1
  img<-png::readPNG(paste("/home/labosea1/image/",sp,".png",sep=""))
  g <- grid::rasterGrob(img, interpolate=TRUE)
  
  data_fake=data.frame(X=seq(0,10,by=1),
                       Y=seq(0,10,by=1))
  p<-ggplot(data_fake,aes(x=X,y=Y)) +
    theme_void() +
    annotation_custom(g) +
    ylim(c(25,100))
  #p
  print(ggpubr::ggarrange(p,p1,nrow=2,heights = c(0.25,0.75)))
  pdf(paste("/home/labosea1/ARTICLE/mitoch/figures/tree_",sp,".pdf",sep=""),width =7.5,height=5)
  print(ggpubr::ggarrange(p,p1,nrow=2,heights = c(0.25,0.75)))
  dev.off()

}

pdf(paste("/home/labosea1/ARTICLE/mitoch/figures/tree.pdf",sep=""),width =20,height=30)
print(
  ggpubr::ggarrange(myplots[[1]],
                    NULL,
                    myplots[[2]],
                    NULL,
                    myplots[[3]],
                    myplots[[4]],
                    NULL,
                    myplots[[5]],
                    NULL,
                    myplots[[6]],
                    myplots[[7]],
                    NULL,
                    myplots[[8]],
                    NULL,
                    myplots[[9]],
                    myplots[[10]],
                    NULL,
                    myplots[[11]],
                    NULL,
                    myplots[[12]],
                    myplots[[13]],
                    NULL,
                    myplots[[14]],
                    NULL,
                    myplots[[15]],
                    myplots[[16]],
                    NULL,
                    myplots[[17]],
                    NULL,
                    myplots[[18]],
                    myplots[[19]],
                    NULL,
                    myplots[[20]],
                    common.legend = T,ncol=5,nrow=7,
                    widths=c(0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3))
)
dev.off()

myplots<-vector('list',20)
a=0
for (sp in SPECIES){
  
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
  a=a+1
  myplots[[a]] <- local({
    p1<-ggtree(tree,col="black",size=1,alpha=0.5) %<+% image_info + 
      #geom_tiplab(hjust=-0.15,size=3,fontface="bold",aes(color=loc)) +
      geom_tippoint(aes(color=loc),size=6,alpha=0.75) +
      theme_tree2() +
      scale_color_manual(name="Location",
                         values=RColorBrewer::brewer.pal(n=4,name="RdBu"))+
      scale_fill_manual(name="Species",values=wesanderson::wes_palette("Zissou1",n=length(levels(factor(image_info$sp))),type="continuous")) +
      xlim(c(0,0.11)) +
      theme(legend.position = "bottom",
            plot.title=element_text(hjust=0.5)) 
    #labs(title="iqtree")
    p1
    img<-png::readPNG(paste("/home/labosea1/image/",sp,".png",sep=""))
    g <- grid::rasterGrob(img, interpolate=TRUE)
    
    data_fake=data.frame(X=seq(0,10,by=1),
                         Y=seq(0,10,by=1))
    p<-ggplot(data_fake,aes(x=X,y=Y)) +
      theme_void() +
      annotation_custom(g) +
      ylim(c(25,100))
    #p
    print(ggpubr::ggarrange(p,p1,nrow=2,heights = c(0.1,0.95)))
  })
  
  
}

pdf(paste("/home/labosea1/ARTICLE/mitoch/figures/tree_scale.pdf",sep=""),width =20,height=30)
print(
  ggpubr::ggarrange(NULL,
                    NULL,
                    myplots[[2]],
                    NULL,
                    myplots[[3]],
                    myplots[[4]],
                    NULL,
                    myplots[[5]],
                    NULL,
                    myplots[[6]],
                    myplots[[7]],
                    NULL,
                    myplots[[8]],
                    NULL,
                    myplots[[9]],
                    myplots[[10]],
                    NULL,
                    myplots[[11]],
                    NULL,
                    myplots[[12]],
                    myplots[[13]],
                    NULL,
                    myplots[[14]],
                    NULL,
                    myplots[[15]],
                    myplots[[16]],
                    NULL,
                    myplots[[17]],
                    NULL,
                    myplots[[18]],
                    myplots[[19]],
                    NULL,
                    myplots[[20]],
                    common.legend = T,ncol=5,nrow=7,
                    widths=c(0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3,0.05,0.3,
                             0.3,0.05,0.3))
)
dev.off()
###

ggplot(data_mitoch[data_mitoch$LENGTH<=10000,],aes(x=GENE,y=TAJIMAD_ALL)) +
  geom_bar(stat='identity',alpha=0.5) +
  facet_wrap(SP~.,nrow=3,scales="free") +
  theme_classic() +
  geom_hline(yintercept=0,col="red",lty=2)

ggplot(data_mitoch[data_mitoch$LENGTH<=10000,],aes(x=LENGTH,y=PI_ALL)) +
  geom_point(aes(col=type_gene),alpha=0.5) +
  facet_wrap(SP~.,nrow=3,scales="free") +
  theme_classic() +
  geom_smooth(method="lm")

data_mitoch[,c(1,2,20:25)] %>%
  filter(GENE != "whole_mitochondrial") %>% 
  reshape2::melt(id.vars=c("SP","GENE")) %>%
  ggplot(aes(x=GENE,y=value)) +
  geom_point(position="jitter",stat='identity',aes(col=variable),alpha=0.4) +
  geom_bar(position="dodge",stat='identity',aes(fill=variable),alpha=0.4) +
  scale_fill_viridis_d(begin=0.1,end=0.9) +
  theme_classic() +
  facet_wrap(SP ~. ,nrow=5,scales="free") +
  ylab("Hudson's Fst") +
  xlab("Species") +
  scale_color_viridis_d(begin=0.1,end=0.9) 


## Ratio


  
####
  
  SPECIES=c(
    "Aboye",
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
    "Styph"
  )
  
  fold_4=c("ct","gt","tc","cc","ac","gc","cg","gg")
  
  
  data_4fold=data.frame(SP=c(NA),
                        GENE=c(NA),
                        count_4fold=c(NA))
  
  for (sp in SPECIES){
    
    setwd(paste("/home/labosea1/MitoZ/ANALYSIS/",sp,"/GENE",sep=""))
    
    for (i in list.files()[grep('_align.fa',list.files())]){
      
      a<-seqinr::read.fasta(i)
      a<-as.character(a[[1]])
      count_4fold=0
      for (j in seq(3,length(a),by=3)){
        test<-paste(a[c(j-2,j-1)],collapse="")
        
        if (test %in% fold_4){
          count_4fold=count_4fold+1
        }
      }
      
      data_4fold=rbind(data_4fold,
                       c(sp,strsplit(i,"_")[[1]][2],count_4fold))
      
      
      
    }
    
  }

  data_4fold=data_4fold[-1,]  

  data_mitoch_4fold=merge(data_mitoch,data_4fold,on=c("SP","GENE"))
data_mitoch_4fold$count_4fold=as.numeric(data_mitoch_4fold$count_4fold)
  
summary(lm(PI_ALL~count_4fold,data=data_mitoch_4fold[data_mitoch_4fold$type_gene=="CDS" & data_mitoch_4fold=="Cjuli",]))

  summary(lm(DXY_LI_GA~count_4fold,data=data_mitoch_4fold[data_mitoch_4fold$type_gene=="CDS" & data_mitoch_4fold$SP=="Dlabr",]))
  plot(data_mitoch_4fold[data_mitoch_4fold$type_gene=="CDS",]$count_4fold,data_mitoch_4fold[data_mitoch_4fold$type_gene=="CDS",]$PI_ALL)  
  
  summary(lm(DXY_LI_GA~count_4fold,data=data_mitoch_4fold[data_mitoch_4fold$type_gene=="CDS",]))
  plot(data_mitoch_4fold[data_mitoch_4fold$type_gene=="CDS",]$count_4fold,data_mitoch_4fold[data_mitoch_4fold$type_gene=="CDS",]$DXY_LI_GA)  
  