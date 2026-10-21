
library('gridExtra')
library('ggplot2')

# setting working directory 
setwd("/home/gospozha/haifa/hiba/sp_align/snp")

VCFsnps <- read.csv('GVCFall_SNPs.table', header = T, na.strings=c("","NA"), sep = "\t") 
VCFindel <- read.csv('GVCFall_INDELs.table', header = T, na.strings=c("","NA"), sep = "\t")
dim(VCFsnps)
#1433670      10
dim(VCFindel)
# 96273     10
VCF <- rbind(VCFsnps, VCFindel)
VCF$Variant <- factor(c(rep("SNPs", dim(VCFsnps)[1]), rep("Indels", dim(VCFindel)[1])))

snps <- '#A9E2E4'
indels <- '#F4CCCA'

DP <- ggplot(VCF, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3) + 
  geom_vline(xintercept=c(10,6200)) 

QD <- ggplot(VCF, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=2, size=0.7)

FS <- ggplot(VCF, aes(x=FS, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(60, 200), size=0.7) + ylim(0,0.1)

MQ <- ggplot(VCF, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=30, size=0.7)

MQRankSum <- ggplot(VCF, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3)

SOR <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(4, 10), size=1, colour = c(snps,indels))

ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) +
  geom_vline(xintercept=c(-10,10,-20,20), size=1, colour = c(snps,snps,indels,indels)) + xlim(-30, 30)

pdf("Diag_plots.pdf", height=20, width=15)
theme_set(theme_gray(base_size = 18))
grid.arrange(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow=4)
dev.off()
