rm(list = ls())
options(stringsAsFactors = F)
library(pheatmap)
library(GSVA)
exp = read.csv("AAV-CBLN1-OE-TPM-oder.csv",header=T,row.names = 1)
AR_gene=c("Ar","Tmprss2","Nkx3-1","Plpp1","Pmepa1","Aldh1a3","Steap4","Fkbp5")
AR_exp = log2(exp[AR_gene,]+1)

exprset=log2(exp+1)
data_tpm=exprset
data_zscore=scale(exprset)
data_zscore <- as.matrix(data_zscore)
data_tpm <- as.matrix(data_tpm)

#if input is counts,then: "kcdf="Possion",others uses default
AR_gene=as.data.frame(AR_gene)
gsvaParam<-gsvaParam(data_zscore, geneSets = AR_gene, minSize = 1)##新版GSVA需要建立gsvapara参数
my_gsva<-gsva(gsvaParam)#GSVA
my_gsva <- data.frame(my_gsva)

AR_exp=rbind(my_gsva,AR_exp)
AR_exp=AR_exp[,c(4:6,1:3)]
colnames(AR_exp)[1:6]=c("AAV_CTL_1","AAV_CTL_2","AAV_CTL_3","AAV_TH_Cbln1_1","AAV_TH_Cbln1_2","AAV_TH_Cbln1_3")

ann_col = data.frame(Sample=factor(c(rep("CTL",3),rep("CBLN1_OE",3))))#创建分组列
row.names(ann_col) = colnames(AR_exp)
ann_color = list(Sample = c(CBLN1_OE="red", CTL="blue")) #定义分组颜色
pheatmap(AR_exp,
         cluster_rows = F,
         cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = T,
         #border_color = NA,
         scale = "row",
         show_rownames =T,
         annotation_col = ann_col,
         annotation_colors=ann_color,
         gaps_row = c(1),
         annotation_legend = TRUE, #表示是否显示注释的图例信息
         
         annotation_names_row = TRUE, annotation_names_col = TRUE)

write.csv(AR_exp,file = "AR_score_heatmap.csv")



