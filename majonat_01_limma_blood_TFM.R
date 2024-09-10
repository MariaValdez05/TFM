rm(list = ls())

library(limma)
library(pheatmap)
library(edgeR)

counts=read.table("COUNTS_genes_MAJONAT_01.tsv",h=T,row.names=1)
info=read.table("info_MAJONAT_01_with_sex_Blood_SINOUTLIERS",h=T, row.names=1)

counts=counts[,names(counts) %in% row.names(info)]
info=info[rownames(info) %in% names(counts),]
info=info[names(counts),]
counts=counts[,rownames(info)]

group=factor(info$GROUP)
tissue=factor(info$TISSUE)
sex=factor(info$SEX)

y=DGEList(counts=counts)
isexpr <- rowSums(cpm(y) > 1) >= 2
y=y[isexpr,keep.lib.size=FALSE]
y=calcNormFactors(y)
dim(y)

combi=combn(unique(c(1,2,3,4)), 2)


mod <- model.matrix(~sex+group)
colnames(mod)=gsub("group","",colnames(mod))

contr.matrix <- makeContrasts(
	

Blood_48hpi_Resistant_vs_Blood_control	=-Blood_control,
Blood_48hpi_Susceptible_vs_Blood_control	=Blood_48hpi_Susceptible-Blood_control,
Blood_48hpi_Resistant_vs_Blood_48hpi_Susceptible=	-Blood_48hpi_Susceptible,

	
levels=colnames(mod))

v=voom(y,mod)
fit=lmFit(v,mod)
fit=contrasts.fit(fit, contrasts=contr.matrix)
fit2=eBayes(fit)
summary(decideTests(fit2))
topTable(fit2)
for (i in  colnames(fit2$coefficients)){
top=topTable(fit2,coef=i,sort="p", n=Inf)
write.table(top,paste(i,"_limma_voom_MAJONAT_01_Blood_corrected_by_sex.txt",sep=""),quote=F)
}

for (i in  colnames(fit2$coefficients)){
top=topTable(fit2,coef=i,sort="p", n=Inf)
genes=rownames(top[which(top$adj.P.Val<0.05),])
term1=strsplit(i,split="_vs_")[[1]][1]
term2=strsplit(i,split="_vs_")[[1]][2]
samples=rownames(subset(info,GROUP==term1 | GROUP==term2))
expr=v$E[genes,samples]
rownames(expr)=do.call(rbind, strsplit(genes, ','))[,2]
if (length(genes) >1) {
#pdf(paste("pheatmap_DE_genes_MAJONAT_01_Lung_corrected_by_sex",i,".pdf",sep=""))
pheatmap(expr,scale="row",annotation_col=info[,c("GROUP","TISSUE","SEX")], border_color = "NA",show_rownames = F, cutree_cols = 2)
#pheatmap(expr,scale="row", border_color = "NA",show_rownames = F, cutree_cols = 2)
#dev.off()
}}

write.table(v$E,"logcpm_MAJONAT_01_Lung.txt",quote=F)
