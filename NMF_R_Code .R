
EAS_count_matrix<-read.table("symbol.txt",header=T,row.names=1)

#convert FPKM to TPM
fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
CGGA_tpm <- apply(EAS_count_matrix,2,fpkmToTpm)
head(CGGA_tpm)  

#log transformation
CGGA_tpm_log <- log2(CGGA_tpm + 1)
rownames(CGGA_tpm_log) <- rownames(CGGA_tpm)
colnames(CGGA_tpm_log) <- colnames(CGGA_tpm)
CGGA_tpm_log <- as.data.frame(CGGA_tpm_log)
head(CGGA_tpm_log)  

#MAD computation and top 3000 genes selection
num=apply(CGGA_tpm_log,1,function(x){length(which(x>0))})
CGGA_tpm_log <- CGGA_tpm_log[which(num>10), ]##至少10个样本表达不为0
EAS3000_countsNorm_log <- as.data.frame(CGGA_tpm_log[order(apply(CGGA_tpm_log,1,mad), decreasing = T)[1:3000],])
EAS3000_countsNorm_log[1:5,1:6]
write.table(EAS3000_countsNorm_log,"MAD3000.txt",sep="\t",quote=FALSE)##for further submap analysis

#NMF analysis
library(NMF)
res2_7 <- nmf(EAS3000_countsNorm_log, 2:7, nrun=30, .options='v')
consensusmap(res2_7, labCol = NA, labRow = NA)
pdf(file="NMF1.pdf",width=12,height=7)
consensusmap(res2_7, labCol = NA, labRow = NA)
dev.off()
pdf(file="NMF2.pdf",width=12,height=10)
plot(res2_7)
dev.off()

coph <- res2_7$measures$cophenetic
coph_diff <- NULL
for (i in 2:length(coph)) {  coph_diff <- c(coph_diff, coph[i-1]-coph[i])}
k.best <- which.max(coph_diff)+1
k.best

#k.beat
est <- nmf(EAS3000_countsNorm_log,3,nrun=10,seed=123456) 

W_mat=basis(est)
W_mat[1:5,]

H_mat=coef(est)
H_mat[,1:5]

group <- as.data.frame(predict(est))
table(group)

write.table(group,file="EAS_NMF_group.txt",sep="\t",quote=FALSE)
write.table(W_mat,file="W_matrix.txt",sep="\t",quote=FALSE)
write.table(H_mat,file="H_matrix.txt",sep="\t",quote=FALSE)

layout(cbind(1, 2))
pdf(file="NMF3.pdf",width=12,height=10)
basismap(est)
dev.off()
pdf(file="NMF4.pdf",width=12,height=10)
coefmap(est)
dev.off()

#PCA
data.pca <- prcomp(EAS3000_countsNorm_log, center=F, scale=F)
EAS_PCA <- as.data.frame(data.pca$rotation)
colnames(group) <- "group"
EAS_PCA <- cbind(group,EAS_PCA)
EAS_PCA$group <- as.factor(EAS_PCA$group)
library(egg)
p1 <- ggplot(EAS_PCA, aes(x = PC1, y = PC2, shape = group, color = group)) + geom_point(size=2)  
p2 <- ggplot(EAS_PCA, aes(x = PC1, y = PC3, shape = group, color = group)) + geom_point(size=2)  
p3 <- ggplot(EAS_PCA, aes(x = PC2, y = PC3, shape = group, color = group)) + geom_point(size=2)  

pdf(file="NMF5.pdf",width=12,height=3)
ggarrange(p1, p2, p3, nrow=1)
dev.off()















