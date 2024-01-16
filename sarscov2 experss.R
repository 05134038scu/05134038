
library(BiocManager)
library(DESeq2)
library(biomaRt)
library(VennDiagram)
library(UpSetR)
#for excrat common data from vennlist
library(ComplexHeatmap)
setwd("/home/yuanpeng/yuan/R sarsCoV2/result")

fr <- read.table(file = "Host.tab", header = TRUE, sep = "\t")

# change remove rowname
myrownames=fr[,1]
fr2=fr
fr2=fr2[,-c(1,2)]
rownames(fr2)=myrownames

# change col name
rsmID <- c("0_r1_notinf","0_r2_notinf","0_r3_notinf","3_r1_notinf","3_r2_notinf","3_r3_notinf","6_r1_notinf","6_r2_notinf","6_r3_notinf","12_r1_notinf","12_r2_notinf","12_r3_notinf","24_r1_notinf","24_r2_notinf","24_r3_notinf","3_r1_inf","3_r2_inf","3_r3_inf","6_r1_inf","6_r2_inf","6_r3_inf","12_r1_inf","12_r2_inf","12_r3_inf","24_r1_inf","24_r2_inf","24_r3_inf")
colnames(fr2)= rsmID

#
fr2=as.matrix(fr2)
h=dim(fr2)
h[1]
num.gene=h[1]

# remove row sum with little expression
# m= row number
m= 0
for(i in 1:num.gene){
  m[i]=mean(fr2[i,])}
# remove low expression
fr2=fr2[m>5,]
# mean made histogram
hist(log(m+1),100)


# check overexpression will showed outside 
plot(log(fr2[,1]+1), log(fr2[,2]+1), pch=".")


#Create DESeq2Dataset object
# first make condition
# install deseq2
condition=c("controls","controls","controls","controls","controls","controls","controls","controls","controls","controls","controls","controls","controls","controls","controls","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected","infected")
dds <- DESeqDataSetFromMatrix(fr2, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)
dds2=dds

#The normalization 
Mnormalized <- counts(dds, normalized=TRUE)
Mnormalized = log(Mnormalized +1)

# intensity filter
mean1=0
for(i in 1:length(Mnormalized[,1])){
  mean1[i]=median(Mnormalized[i,])}
Mnormalized2=Mnormalized[mean1>quantile(mean1,0.3),]

# variance filter
var1=0
for(i in 1:length(Mnormalized2[,1])){
  var1[i]=sd(Mnormalized2[i,])}
Mnormalized3=Mnormalized2[var1>quantile(var1,0.2),]

#remove mock 3h-24h  cluster前再移除
Mnormalized3= Mnormalized3[,!colnames(Mnormalized3) %in% c("3_r1_notinf","3_r2_notinf","3_r3_notinf","6_r1_notinf","6_r2_notinf","6_r3_notinf","12_r1_notinf","12_r2_notinf","12_r3_notinf","24_r1_notinf","24_r2_notinf","24_r3_notinf")]


# clustering 
method <- "average"
distance <- "correlation"
clustering <- pvclust((Mnormalized3),method.hclust="average", method.dist="correlation")
plot(clustering)


#remove some outlier first
Mnormalized3= Mnormalized3[,!colnames(Mnormalized3) %in% c("0_r1_notinf", "3_r2_inf", "12_r2_inf", "3_r1_inf")]
#cluster again
clustering <- pvclust((Mnormalized3),method.hclust="average", method.dist="correlation")
plot(clustering)

# example scatterplot of normlaized data
plot(Mnormalized3[,1],Mnormalized3[,2],pch=".")

# remove outlier from raw data so can do deseq
fr3=fr2
fr3= fr3[,!colnames(fr3) %in% c("3_r1_notinf","3_r2_notinf","3_r3_notinf","6_r1_notinf","6_r2_notinf","6_r3_notinf","12_r1_notinf","12_r2_notinf","12_r3_notinf","24_r1_notinf","24_r2_notinf","24_r3_notinf","0_r1_notinf", "3_r2_inf", "12_r2_inf", "3_r1_inf")]
fr3=as.matrix(fr3)

#>>>>>>>>>>>>.new total>>>>>>>>>>>>
m= 0
for(i in 1:nrow(fr3)){
  m[i]=sum(fr3[i,])}
# remove low expression
fr3=fr3[m>11,]
fr4=as.data.frame(fr3+1) 
write.csv(fr4, file = "mrna_fr3.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
#>>>>>>>>>>>>.new total>>>>>>>>>>>>


# 6 hrs vs control
cfr3=fr3[,c(1,2,3,4,5,6)]
condition=c("controls","controls","infected","infected","infected","infected")

#Create DESeq2Dataset object

dds <- DESeqDataSetFromMatrix(cfr3, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)
res <- results(dds)

str(res)
res$padj[is.na(res$padj)]=1

# select up down regulate
upregulated6h=c()
downregulated6h=c()
rw=rownames(fr3)
for(i in 1: length(res$pvalue)){
  if(res$padj[i]<=0.001){
    if(res$log2FoldChange[i] > 0.7)
    {upregulated6h=rbind(upregulated6h, c(rw[i],res$log2FoldChange[i],res$padj[i]))}
    if(res$log2FoldChange[i] < -0.7)
    {downregulated6h=rbind(downregulated6h, c(rw[i],res$log2FoldChange[i],res$padj[i]))}
  }}

upregulated6h=data.frame(upregulated6h)
colnames(upregulated6h)=c("ID", "log2FoldChange", "padj")
downregulated6h= data.frame(downregulated6h)
colnames(downregulated6h)=c("ID", "log2FoldChange", "padj")
upregulated6h$log2FoldChange=as.numeric(upregulated6h$log2FoldChange)
upregulated6h$padj=as.numeric(upregulated6h$padj)
downregulated6h$log2FoldChange=as.numeric(downregulated6h$log2FoldChange)
downregulated6h$padj=as.numeric(downregulated6h$padj)



#volcano plot total
temp_res= cbind(rw,res$log2FoldChange,res$padj)
colnames(temp_res)=c("ID", "log2FoldChange", "padj")
temp_res=as.data.frame(temp_res)
temp_res$log2FoldChange= as.numeric(temp_res$log2FoldChange)
temp_res$padj= as.numeric(temp_res$padj)

temp_res$diffexpressed= "No"
temp_res$diffexpressed[temp_res$log2FoldChange>0.05 & temp_res$padj<0.05]="Up"
temp_res$diffexpressed[temp_res$log2FoldChange<0.05 & temp_res$padj<0.05]="Down"

temp_res_vol= ggplot(temp_res, aes(x=temp_res$log2FoldChange, y=-log10(temp_res$padj)))+ geom_point() + theme_minimal()
temp_res_vol2= temp_res_vol+ geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")
temp_res_vol3 = temp_res_vol2 + scale_color_manual(values=mycolors)

ggplot(temp_res, aes(x=temp_res_vol3$data$log2FoldChange, y=-log10(temp_res_vol3$data$padj), col=temp_res_vol3$data$diffexpressed))+ 
  geom_point() + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")+ 
  theme_minimal()+labs(title= "Volcano_plot_6h", x = "log2FoldChange", y = "-log10(padj)", colour="Diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#volcano plot
temp_6_updw=rbind(upregulated6h,downregulated6h)
temp_6_updw_vol= ggplot(temp_6_updw, aes(x=temp_6_updw$log2FoldChange, y=-log10(temp_6_updw$padj)))+ geom_point() + theme_minimal()
temp_6_updw_vol2 = temp_6_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

temp_6_updw$diffexpressed= "NO"
temp_6_updw$diffexpressed[temp_6_updw$log2FoldChange>0.05 & temp_6_updw$padj<0.05]="up"
temp_6_updw$diffexpressed[temp_6_updw$log2FoldChange<0.05 & temp_6_updw$padj<0.05]="down"

temp_6_updw_vol= ggplot(temp_6_updw, aes(x=temp_6_updw$log2FoldChange, y=-log10(temp_6_updw$padj), col=temp_6_updw$diffexpressed))+ geom_point() + theme_minimal()
temp_6_updw_vol2 = temp_6_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_6_updw_vol3=temp_6_updw_vol2+ scale_color_manual(values=mycolors)

ggplot(temp_6_updw, aes(x=temp_6_updw_vol3$data$log2FoldChange, y=-log10(temp_6_updw_vol3$data$padj), col=temp_6_updw_vol3$data$diffexpressed))+ 
  geom_point() + theme_minimal()+labs(title= "volcano_plot_6h", x = "log2FoldChange", y = "-log10(padj)", colour="diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# 12 hrs vs control
cfr3=fr3[,c(1,2,7,8)]
condition=c("controls","controls","infected","infected")

#Create 12 hr DESeq2Dataset object

dds <- DESeqDataSetFromMatrix(cfr3, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)
res <- results(dds)

str(res)
res$padj[is.na(res$padj)]=1

# select 12hr up down regulate
upregulated12h=c()
downregulated12h=c()
rw=rownames(fr3)
for(i in 1: length(res$pvalue)){
   if(res$padj[i]<=0.001){
     if(res$log2FoldChange[i] > 0.7)
     {upregulated12h=rbind(upregulated12h, c(rw[i],res$log2FoldChange[i],res$padj[i]))}
     if(res$log2FoldChange[i] < -0.7)
     {downregulated12h=rbind(downregulated12h, c(rw[i],res$log2FoldChange[i],res$padj[i]))}
}}
 
upregulated12h=data.frame(upregulated12h)
colnames(upregulated12h)=c("ID", "log2FoldChange", "padj")
downregulated12h= data.frame(downregulated12h)
colnames(downregulated12h)=c("ID", "log2FoldChange", "padj")
upregulated12h$log2FoldChange=as.numeric(upregulated12h$log2FoldChange)
upregulated12h$padj=as.numeric(upregulated12h$padj)
downregulated12h$log2FoldChange=as.numeric(downregulated12h$log2FoldChange)
downregulated12h$padj=as.numeric(downregulated12h$padj)

#volcano plot total
temp_res= cbind(rw,res$log2FoldChange,res$padj)
colnames(temp_res)=c("ID", "log2FoldChange", "padj")
temp_res=as.data.frame(temp_res)
temp_res$log2FoldChange= as.numeric(temp_res$log2FoldChange)
temp_res$padj= as.numeric(temp_res$padj)
 
temp_res$diffexpressed= "No"
temp_res$diffexpressed[temp_res$log2FoldChange>0.05 & temp_res$padj<0.05]="Up"
temp_res$diffexpressed[temp_res$log2FoldChange<0.05 & temp_res$padj<0.05]="Down"


temp_res_vol= ggplot(temp_res, aes(x=temp_res$log2FoldChange, y=-log10(temp_res$padj)))+ geom_point() + theme_minimal()
temp_res_vol2= temp_res_vol+ geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")
temp_res_vol3 = temp_res_vol2 + scale_color_manual(values=mycolors)

ggplot(temp_res, aes(x=temp_res_vol3$data$log2FoldChange, y=-log10(temp_res_vol3$data$padj), col=temp_res_vol3$data$diffexpressed))+ 
  geom_point() + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")+ 
  theme_minimal()+labs(title= "Volcano_plot_12h", x = "log2FoldChange", y = "-log10(padj)", colour="Diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))


#>>>>>>>>>>>>>>>>
#volcano plot 
temp_12_updw=rbind(upregulated12h,downregulated12h)
temp_12_updw_vol= ggplot(temp_12_updw, aes(x=temp_12_updw$log2FoldChange, y=-log10(temp_12_updw$padj)))+ geom_point() + theme_minimal()
temp_12_updw_vol2 = temp_12_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

temp_12_updw$diffexpressed= "NO"
temp_12_updw$diffexpressed[temp_12_updw$log2FoldChange>0.05 & temp_12_updw$padj<0.05]="up"
temp_12_updw$diffexpressed[temp_12_updw$log2FoldChange<0.05 & temp_12_updw$padj<0.05]="down"

temp_12_updw_vol= ggplot(temp_12_updw, aes(x=temp_12_updw$log2FoldChange, y=-log10(temp_12_updw$padj), col=temp_12_updw$diffexpressed))+ geom_point() + theme_minimal()
temp_12_updw_vol2 = temp_12_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_12_updw_vol3=temp_12_updw_vol2+ scale_color_manual(values=mycolors)

ggplot(temp_12_updw, aes(x=temp_12_updw_vol3$data$log2FoldChange, y=-log10(temp_12_updw_vol3$data$padj), col=temp_12_updw_vol3$data$diffexpressed))+ 
  geom_point() + theme_minimal()+labs(title= "volcano_plot_12h", x = "log2FoldChange", y = "-log10(padj)", colour="diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))
#>>>>>>>>>>>>>>>>>>>>>>>


# 24 hrs vs control
cfr3=fr3[,c(1,2,9,10,11)]
condition=c("controls","controls","infected","infected","infected")

#Create 24hr DESeq2Dataset object

dds <- DESeqDataSetFromMatrix(cfr3, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)
res <- results(dds)

str(res)
res$padj[is.na(res$padj)]=1

# select 24hr up down regulate
upregulated24h=c()
downregulated24h=c()
rw=rownames(fr3)
for(i in 1: length(res$pvalue)){
  if(res$padj[i]<=0.001){
    if(res$log2FoldChange[i] > 0.7)
    {upregulated24h=rbind(upregulated24h, c(rw[i],res$log2FoldChange[i],res$padj[i]))}
    if(res$log2FoldChange[i] < -0.7)
    {downregulated24h=rbind(downregulated24h, c(rw[i],res$log2FoldChange[i],res$padj[i]))}
    }}

upregulated24h=data.frame(upregulated24h)
colnames(upregulated24h)=c("ID", "log2FoldChange", "padj")
downregulated24h= data.frame(downregulated24h)
colnames(downregulated24h)=c("ID", "log2FoldChange", "padj")
upregulated24h$log2FoldChange=as.numeric(upregulated24h$log2FoldChange)
upregulated24h$padj=as.numeric(upregulated24h$padj)
downregulated24h$log2FoldChange=as.numeric(downregulated24h$log2FoldChange)
downregulated24h$padj=as.numeric(downregulated24h$padj)

#volcano plot total
temp_res= cbind(rw,res$log2FoldChange,res$padj)
colnames(temp_res)=c("ID", "log2FoldChange", "padj")
temp_res=as.data.frame(temp_res)
temp_res$log2FoldChange= as.numeric(temp_res$log2FoldChange)
temp_res$padj= as.numeric(temp_res$padj)

temp_res$diffexpressed= "No"
temp_res$diffexpressed[temp_res$log2FoldChange>0.05 & temp_res$padj<0.05]="Up"
temp_res$diffexpressed[temp_res$log2FoldChange<0.05 & temp_res$padj<0.05]="Down"


temp_res_vol= ggplot(temp_res, aes(x=temp_res$log2FoldChange, y=-log10(temp_res$padj)))+ geom_point() + theme_minimal()
temp_res_vol2= temp_res_vol+ geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")
temp_res_vol3 = temp_res_vol2 + scale_color_manual(values=mycolors)

ggplot(temp_res, aes(x=temp_res_vol3$data$log2FoldChange, y=-log10(temp_res_vol3$data$padj), col=temp_res_vol3$data$diffexpressed))+ 
  geom_point() + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")+ 
  theme_minimal()+labs(title= "Volcano_plot_24h", x = "log2FoldChange", y = "-log10(padj)", colour="Diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#volcano plot without no differential express
temp_24_updw=rbind(upregulated24h,downregulated24h)
temp_24_updw_vol= ggplot(temp_24_updw, aes(x=temp_24_updw$log2FoldChange, y=-log10(temp_24_updw$padj)))+ geom_point() + theme_minimal()
temp_24_updw_vol2 = temp_24_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

temp_24_updw$diffexpressed= "No"
temp_24_updw$diffexpressed[temp_24_updw$log2FoldChange>0.05 & temp_24_updw$padj<0.05]="Up"
temp_24_updw$diffexpressed[temp_24_updw$log2FoldChange<0.05 & temp_24_updw$padj<0.05]="Down"

temp_24_updw_vol= ggplot(temp_24_updw, aes(x=temp_24_updw$log2FoldChange, y=-log10(temp_24_updw$padj), col=temp_24_updw$diffexpressed))+ geom_point() + theme_minimal()
temp_24_updw_vol2 = temp_24_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

names(mycolors) <- c("Down", "Up", "No")
temp_24_updw_vol3=temp_24_updw_vol2+ scale_color_manual(values=mycolors)

ggplot(temp_24_updw, aes(x=temp_24_updw_vol3$data$log2FoldChange, y=-log10(temp_24_updw_vol3$data$padj), col=temp_24_updw_vol3$data$diffexpressed))+ 
  geom_point() + theme_minimal()+labs(title= "volcano_plot_24h", x = "log2FoldChange", y = "-log10(padj)", colour="diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))

# 0920 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
upgene= bind_rows(upregulated6h,upregulated12h,upregulated24h)
upgene2 <- upgene %>%
  distinct(ID)
dwgene= bind_rows(downregulated6h,downregulated12h,downregulated24h)
dwgene2 <- dwgene %>%
  distinct(ID)

upgn_fr=fr3[row.names(fr3) %in% upgene2$ID, ]  
dwgn_fr=fr3[row.names(fr3) %in% dwgene2$ID, ]  
updwgn_fr=rbind(upgn_fr,dwgn_fr)


write.csv(upgn_fr, file = "upgn_fr_mrna.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
write.csv(dwgn_fr, file = "dwgn_fr_mrna.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
write.csv(updwgn_fr, file = "updwgn_fr_mrna.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)

# 轉乘gn id 
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_ID=biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name'),
                       filters = 'ensembl_gene_id', 
                       values = rownames(updwgn_fr),
                       mart = mart)
#查爲什饃有些ENST不見
miss= updwgn_fr[!rownames(updwgn_fr) %in% gene_ID$ensembl_gene_id,]
#從 cir那邊的方法用的
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
updwgn_fr <- data.frame("rowname" = rownames(updwgn_fr), updwgn_fr)
updwgn_fr2= merge(updwgn_fr,gene_ID, by.x = 'rowname', by.y = 'ensembl_gene_id', all = FALSE)
rownames(updwgn_fr2) <- updwgn_fr2$rowname



colnames(updwgn_fr2)= c("ENSG","0_r2_notinf","0_r3_notinf","3_r3_inf","6_r1_inf","6_r2_inf","6_r3_inf",
                        "12_r1_inf","12_r3_inf","24_r1_inf","24_r2_inf","24_r3_inf","hgnc","gn_name" )

write.csv(updwgn_fr2, file = "updwgn_fr_cir2.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)


# 0920 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>



####>>>>>>>>>>>>>>>>>>>>>>>>>>>>.


#connect to ensembl  biomart
# h1 = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", GRCh = 37)

# make object to store list after it convert to gene name
genenames_6h_up = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = upregulated6h, mart = ensembl)
genenames_12h_up = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = upregulated12h, mart = ensembl)
genenames_24h_up = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = upregulated24h, mart = ensembl)

genenames_6h_dw = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = downregulated6h, mart = ensembl)
genenames_12h_dw = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = downregulated12h, mart = ensembl)
genenames_24h_dw = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = downregulated24h, mart = ensembl)

# venn diagram  make list 
venn_uplist= list ("inf_6hr_up"=genenames_6h_up$hgnc_symbol,"inf_12hr_up"=genenames_12h_up$hgnc_symbol,"inf_24hr_up"=genenames_24h_up$hgnc_symbol)
venn.diagram(venn_uplist, filename = 'venn_up.png', imagetype = 'png', fill = c('red', 'blue','green'), alpha = 0.50, cat.col = c('red', 'blue', 'green'), 
             col = c('red', 'blue', 'green'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')

venn_dwlist= list ("inf_6hr_dw"=genenames_6h_dw$hgnc_symbol,"inf_12hr_dw"=genenames_12h_dw$hgnc_symbol,"inf_24hr_dw"=genenames_24h_dw$hgnc_symbol)
venn.diagram(venn_dwlist, filename = 'venn_dw.png', imagetype = 'png', fill = c('red', 'blue','green'), alpha = 0.50, cat.col = c('red', 'blue', 'green'), 
             col = c('red', 'blue', 'green'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')

#upset plot -upregulate 
#venn_uplist3 <- UpSetR::fromList(venn_uplist) 
upset(venn_uplist3, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE))
# 2way to do this
upset(fromList(venn_uplist), mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE))

pdf("upsetplot_up.pdf")
upset(fromList(venn_uplist), mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(TRUE,T))
dev.off()

pdf("upsetplot_up2.pdf")
upset(fromList(venn_uplist), mb.ratio = c(0.5, 0.5))
dev.off()

# upset plot -dwregulate 
pdf("upsetplot_dw.pdf")
upset(fromList(venn_dwlist), mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(TRUE,T))
dev.off()

pdf("upsetplot_dw2.pdf")
upset(fromList(venn_dwlist), mb.ratio = c(0.5, 0.5))
dev.off()

#extract common data
comdata_up= make_comb_mat(venn_uplist)
set_name(comdata_up)
comb_name(comdata_up)
set_size(comdata_up)
comb_size(comdata_up)

# extract common up data
com_upgene_6v12 <- extract_comb(comdata_up, "110")
com_upgene_6v24 = extract_comb(comdata_up, "101")
com_upgene_12v24 = extract_comb(comdata_up, "011")
com_upgene_6v12v24 = extract_comb(comdata_up, "111")
# extract common down data
comdata_dw= make_comb_mat(venn_dwlist)
com_dwgene_6v12 = extract_comb(comdata_dw, "110")
com_dwgene_6v24 = extract_comb(comdata_dw, "101")
com_dwgene_12v24 = extract_comb(comdata_dw, "011")
com_dwgene_6v12v24 = extract_comb(comdata_dw, "111")

#Gene set enrichment analysis
# genelist from up 
gnlist_gnname_6h_up=gprofiler2::gost(genenames_6h_up[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_12h_up=gprofiler2::gost(genenames_12h_up[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_24h_up=gprofiler2::gost(genenames_24h_up[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
# only receive $result from gprofiler
gnlist_gnname_6h_up= gnlist_gnname_6h_up$result
gnlist_gnname_12h_up= gnlist_gnname_12h_up$result
gnlist_gnname_24h_up= gnlist_gnname_24h_up$result

gnlist_com_upgene_6v12= gprofiler2::gost(com_upgene_6v12, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_com_upgene_6v24= gprofiler2::gost(com_upgene_6v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)  
gnlist_com_upgene_12v24= gprofiler2::gost(com_upgene_12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)   
gnlist_com_upgene_6v12v24= gprofiler2::gost(com_upgene_6v12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)

gnlist_com_upgene_6v12= gnlist_com_upgene_6v12$result
gnlist_com_upgene_6v24= gnlist_com_upgene_6v24$result
gnlist_com_upgene_12v24= gnlist_com_upgene_12v24$result
gnlist_com_upgene_6v12v24= gnlist_com_upgene_6v12v24$result

# select only 10-1000 of term size
gnlist_gnname_6h_up2<- data.frame()
for (i in 1:length(gnlist_gnname_6h_up$term_size)){
  if(gnlist_gnname_6h_up$term_size[i]<1000 && gnlist_gnname_6h_up$term_size[i]>10){
    print(i)
    print(gnlist_gnname_6h_up$term_size[i])
    gnlist_gnname_6h_up2=rbind(gnlist_gnname_6h_up2,c(gnlist_gnname_6h_up$term_size[i],gnlist_gnname_6h_up$source[i],
                                                      gnlist_gnname_6h_up$p_value[i],gnlist_gnname_6h_up$term_id[i]
                                                      ,gnlist_gnname_6h_up$term_name[i],gnlist_gnname_6h_up$query_size[i]
                                                      ,gnlist_gnname_6h_up$intersection_size[i],gnlist_gnname_6h_up$intersection[i]))
  }
}

colnames(gnlist_gnname_6h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
# order by p value  order(as.numeric) let 0 e-4 become numeric and can be ordered
gnlist_gnname_6h_up2=gnlist_gnname_6h_up2[order(as.numeric(gnlist_gnname_6h_up2$p_value)),]
# change num size of p value to become only 5 number > 0.00001
gnlist_gnname_6h_up2$p_value=round(as.numeric(gnlist_gnname_6h_up2$p_value) , digits = 5)
write.csv(gnlist_gnname_6h_up2, file = "gnlist_gnname_6h_up2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


gnlist_gnname_12h_up2<- data.frame()
for (i in 1:length(gnlist_gnname_12h_up$term_size)){
  if(gnlist_gnname_12h_up$term_size[i]<1000 && gnlist_gnname_12h_up$term_size[i]>10){
    print(i)
    print(gnlist_gnname_12h_up$term_size[i])
    gnlist_gnname_12h_up2=rbind(gnlist_gnname_12h_up2,c(gnlist_gnname_12h_up$term_size[i],gnlist_gnname_12h_up$source[i],
                                                        gnlist_gnname_12h_up$p_value[i],gnlist_gnname_12h_up$term_id[i],
                                                        gnlist_gnname_12h_up$term_name[i],gnlist_gnname_12h_up$query_size[i]
                                                        ,gnlist_gnname_12h_up$intersection_size[i],gnlist_gnname_12h_up$intersection[i]))
  }
}

colnames(gnlist_gnname_12h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_gnname_12h_up2=gnlist_gnname_12h_up2[order(as.numeric(gnlist_gnname_12h_up2$p_value)),]
gnlist_gnname_12h_up2$p_value=round(as.numeric(gnlist_gnname_12h_up2$p_value) , digits = 5)
write.csv(gnlist_gnname_12h_up2, file = "gnlist_gnname_12h_up2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_gnname_24h_up2<- data.frame()
for (i in 1:length(gnlist_gnname_24h_up$term_size)){
  if(gnlist_gnname_24h_up$term_size[i]<1000 && gnlist_gnname_24h_up$term_size[i]>10){
    print(i)
    print(gnlist_gnname_24h_up$term_size[i])
    gnlist_gnname_24h_up2=rbind(gnlist_gnname_24h_up2,c(gnlist_gnname_24h_up$term_size[i],gnlist_gnname_24h_up$source[i],
                                                        gnlist_gnname_24h_up$p_value[i],gnlist_gnname_24h_up$term_id[i],
                                                        gnlist_gnname_24h_up$term_name[i],gnlist_gnname_24h_up$query_size[i]
                                                        ,gnlist_gnname_24h_up$intersection_size[i],gnlist_gnname_24h_up$intersection[i]))
  }
}

colnames(gnlist_gnname_24h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_gnname_24h_up2=gnlist_gnname_24h_up2[order(as.numeric(gnlist_gnname_24h_up2$p_value)),]
gnlist_gnname_24h_up2$p_value=round(as.numeric(gnlist_gnname_24h_up2$p_value) , digits = 5)
write.csv(gnlist_gnname_24h_up2, file = "gnlist_gnname_24h_up2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_com_6v12h_up2<- data.frame()  
for (i in 1:length(gnlist_com_upgene_6v12$term_size)){
  if(gnlist_com_upgene_6v12$term_size[i]<1000 && gnlist_com_upgene_6v12$term_size[i]>10){
    print(i)
    print(gnlist_com_upgene_6v12$term_size[i])
    gnlist_com_6v12h_up2=rbind(gnlist_com_6v12h_up2,c(gnlist_com_upgene_6v12$term_size[i],gnlist_com_upgene_6v12$source[i],
                                                      gnlist_com_upgene_6v12$p_value[i],gnlist_com_upgene_6v12$term_id[i],
                                                      gnlist_com_upgene_6v12$term_name[i],gnlist_com_upgene_6v12$query_size[i]
                                                      ,gnlist_com_upgene_6v12$intersection_size[i],gnlist_com_upgene_6v12$intersection[i]))
  }
}

colnames(gnlist_com_6v12h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_6v12h_up2=gnlist_com_6v12h_up2[order(as.numeric(gnlist_com_6v12h_up2$p_value)),]
gnlist_com_6v12h_up2$p_value=round(as.numeric(gnlist_com_6v12h_up2$p_value) , digits = 5)
write.csv(gnlist_com_6v12h_up2, file = "gnlist_com_6v12h_up2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_com_6v24h_up2<- data.frame()
for (i in 1:length(gnlist_com_upgene_6v24$term_size)){
  if(gnlist_com_upgene_6v24$term_size[i]<1000 && gnlist_com_upgene_6v24$term_size[i]>10){
    print(i)
    print(gnlist_com_upgene_6v24$term_size[i])
    gnlist_com_6v24h_up2=rbind(gnlist_com_6v24h_up2,c(gnlist_com_upgene_6v24$term_size[i],gnlist_com_upgene_6v24$source[i],
                                                      gnlist_com_upgene_6v24$p_value[i],gnlist_com_upgene_6v24$term_id[i],
                                                      gnlist_com_upgene_6v24$term_name[i],gnlist_com_upgene_6v24$query_size[i]
                                                      ,gnlist_com_upgene_6v24$intersection_size[i],gnlist_com_upgene_6v24$intersection[i]))
  }
}

colnames(gnlist_com_6v24h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_6v24h_up2$p_value=round(as.numeric(gnlist_com_6v24h_up2$p_value) , digits = 5)
gnlist_com_6v24h_up2=gnlist_com_6v24h_up2[order(as.numeric(gnlist_com_6v24h_up2$p_value)),]
write.csv(gnlist_com_6v24h_up2, file = "gnlist_com_6v24h_up2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_com_12v24h_up2<- data.frame()
for (i in 1:length(gnlist_com_upgene_12v24$term_size)){
  if(gnlist_com_upgene_12v24$term_size[i]<1000 && gnlist_com_upgene_12v24$term_size[i]>10){
    print(i)
    print(gnlist_com_upgene_12v24$term_size[i])
    gnlist_com_12v24h_up2=rbind(gnlist_com_12v24h_up2,c(gnlist_com_upgene_12v24$term_size[i],gnlist_com_upgene_12v24$source[i],
                                                        gnlist_com_upgene_12v24$p_value[i],gnlist_com_upgene_12v24$term_id[i],
                                                        gnlist_com_upgene_12v24$term_name[i],gnlist_com_upgene_12v24$query_size[i]
                                                        ,gnlist_com_upgene_12v24$intersection_size[i],gnlist_com_upgene_12v24$intersection[i]))
  }
}

colnames(gnlist_com_12v24h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_12v24h_up2$p_value=round(as.numeric(gnlist_com_12v24h_up2$p_value) , digits = 5)
gnlist_com_12v24h_up2=gnlist_com_12v24h_up2[order(as.numeric(gnlist_com_12v24h_up2$p_value)),]
write.csv(gnlist_com_12v24h_up2, file = "gnlist_com_12v24h_up2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_com_6v12v24h_up2<- data.frame()
for (i in 1:length(gnlist_com_upgene_6v12v24$term_size)){
  if(gnlist_com_upgene_6v12v24$term_size[i]<1000 && gnlist_com_upgene_6v12v24$term_size[i]>10){
    print(i)
    print(gnlist_com_upgene_6v12v24$term_size[i])
    gnlist_com_6v12v24h_up2=rbind(gnlist_com_6v12v24h_up2,c(gnlist_com_upgene_6v12v24$term_size[i],gnlist_com_upgene_6v12v24$source[i],
                                                            gnlist_com_upgene_6v12v24$p_value[i],gnlist_com_upgene_6v12v24$term_id[i],
                                                            gnlist_com_upgene_6v12v24$term_name[i],gnlist_com_upgene_6v12v24$query_size[i]
                                                            ,gnlist_com_upgene_6v12v24$intersection_size[i],gnlist_com_upgene_6v12v24$intersection[i]))
  }
}

colnames(gnlist_com_6v12v24h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_6v12v24h_up2$p_value=round(as.numeric(gnlist_com_6v12v24h_up2$p_value) , digits = 5)
gnlist_com_6v12v24h_up2=gnlist_com_6v12v24h_up2[order(as.numeric(gnlist_com_6v12v24h_up2$p_value)),]
write.csv(gnlist_com_6v12v24h_up2, file = "gnlist_com_6v12v24h_up2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


# old version of gprofiler2
#gnlist_com_upgene_6v12=gprofiler(com_upgene_6v12, correction_method="fdr", min_set_size=20, max_set_size=300,hier_filtering="strong", src_filter="GO")
#gnlist_com_upgene_6v24=gprofiler(com_upgene_6v24, correction_method="fdr", min_set_size=20, max_set_size=300,hier_filtering="strong", src_filter="GO")
#gnlist_com_upgene_12v24=gprofiler(com_upgene_12v24, correction_method="fdr", min_set_size=20, max_set_size=300,hier_filtering="strong", src_filter="GO")
#gnlist_com_upgene_6v12v24=gprofiler(com_upgene_6v12v24, correction_method="fdr", min_set_size=20, max_set_size=300,hier_filtering="strong", src_filter="GO")

# genelist from dw
#gnlist_gnname_6h_dw=gprofiler(genenames_6h_dw, correction_method="fdr", min_set_size=20, max_set_size=300,hier_filtering="strong", src_filter="GO")

gnlist_gnname_6h_dw=gprofiler2::gost(genenames_6h_dw[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_12h_dw=gprofiler2::gost(genenames_12h_dw[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_24h_dw=gprofiler2::gost(genenames_24h_dw[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
# only receive $result from gprofiler
gnlist_gnname_6h_dw= gnlist_gnname_6h_dw$result
gnlist_gnname_12h_dw= gnlist_gnname_12h_dw$result
gnlist_gnname_24h_dw= gnlist_gnname_24h_dw$result

gnlist_com_dwgene_6v12= gprofiler2::gost(com_dwgene_6v12, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_com_dwgene_6v24= gprofiler2::gost(com_dwgene_6v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_com_dwgene_12v24= gprofiler2::gost(com_dwgene_12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE) 
gnlist_com_dwgene_6v12v24= gprofiler2::gost(com_dwgene_6v12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)

gnlist_com_dwgene_6v12= gnlist_com_dwgene_6v12$result
gnlist_com_dwgene_6v24= gnlist_com_dwgene_6v24$result
gnlist_com_dwgene_12v24= gnlist_com_dwgene_12v24$result
gnlist_com_dwgene_6v12v24= gnlist_com_dwgene_6v12v24$result

# select only 10<term size<1000
gnlist_gnname_6h_dw2<- data.frame()
for (i in 1:length(gnlist_gnname_6h_dw$term_size)){
  if(gnlist_gnname_6h_dw$term_size[i]<1000 && gnlist_gnname_6h_dw$term_size[i]>10){
    print(i)
    print(gnlist_gnname_6h_dw$term_size[i])
    gnlist_gnname_6h_dw2=rbind(gnlist_gnname_6h_dw2,c(gnlist_gnname_6h_dw$term_size[i],gnlist_gnname_6h_dw$source[i],
                                                      gnlist_gnname_6h_dw$p_value[i],gnlist_gnname_6h_dw$term_id[i],
                                                      gnlist_gnname_6h_dw$term_name[i],gnlist_gnname_6h_dw$query_size[i]
                                                      ,gnlist_gnname_6h_dw$intersection_size[i],gnlist_gnname_6h_dw$intersection[i]))
  }
}

colnames(gnlist_gnname_6h_dw2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_gnname_6h_dw2$p_value=round(as.numeric(gnlist_gnname_6h_dw2$p_value) , digits = 5)
gnlist_gnname_6h_dw2=gnlist_gnname_6h_dw2[order(as.numeric(gnlist_gnname_6h_dw2$p_value)),]
write.csv(gnlist_gnname_6h_dw2, file = "gnlist_gnname_6h_dw2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_gnname_12h_dw2<- data.frame()
for (i in 1:length(gnlist_gnname_12h_dw$term_size)){
  if(gnlist_gnname_12h_dw$term_size[i]<1000 && gnlist_gnname_12h_dw$term_size[i]>10){
    print(i)
    print(gnlist_gnname_12h_dw$term_size[i])
    gnlist_gnname_12h_dw2=rbind(gnlist_gnname_12h_dw2,c(gnlist_gnname_12h_dw$term_size[i],gnlist_gnname_12h_dw$source[i],
                                                        gnlist_gnname_12h_dw$p_value[i],gnlist_gnname_12h_dw$term_id[i],
                                                        gnlist_gnname_12h_dw$term_name[i],gnlist_gnname_12h_dw$query_size[i]
                                                        ,gnlist_gnname_12h_dw$intersection_size[i],gnlist_gnname_12h_dw$intersection[i]))
  }
}

colnames(gnlist_gnname_12h_dw2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_gnname_12h_dw2$p_value=round(as.numeric(gnlist_gnname_12h_dw2$p_value) , digits = 5)
gnlist_gnname_12h_dw2=gnlist_gnname_12h_dw2[order(as.numeric(gnlist_gnname_12h_dw2$p_value)),]
write.csv(gnlist_gnname_12h_dw2, file = "gnlist_gnname_12h_dw2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_gnname_24h_dw2<- data.frame()
for (i in 1:length(gnlist_gnname_24h_dw$term_size)){
  if(gnlist_gnname_24h_dw$term_size[i]<1000 && gnlist_gnname_24h_dw$term_size[i]>10){
    print(i)
    print(gnlist_gnname_24h_dw$term_size[i])
    gnlist_gnname_24h_dw2=rbind(gnlist_gnname_24h_dw2,c(gnlist_gnname_24h_dw$term_size[i],gnlist_gnname_24h_dw$source[i],
                                                        gnlist_gnname_24h_dw$p_value[i],gnlist_gnname_24h_dw$term_id[i],
                                                        gnlist_gnname_24h_dw$term_name[i],gnlist_gnname_24h_dw$query_size[i]
                                                        ,gnlist_gnname_24h_dw$intersection_size[i],gnlist_gnname_24h_dw$intersection[i]))
  }
}

colnames(gnlist_gnname_24h_dw2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_gnname_24h_dw2$p_value=round(as.numeric(gnlist_gnname_24h_dw2$p_value) , digits = 5)
gnlist_gnname_24h_dw2=gnlist_gnname_24h_dw2[order(as.numeric(gnlist_gnname_24h_dw2$p_value)),]
write.csv(gnlist_gnname_24h_dw2, file = "gnlist_gnname_24h_dw2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_com_6v12h_dw2<- data.frame()
for (i in 1:length(gnlist_com_dwgene_6v12$term_size)){
  if(gnlist_com_dwgene_6v12$term_size[i]<1000 && gnlist_com_dwgene_6v12$term_size[i]>10){
    print(i)
    print(gnlist_com_dwgene_6v12$term_size[i])
    gnlist_com_6v12h_dw2=rbind(gnlist_com_6v12h_dw2,c(gnlist_com_dwgene_6v12$term_size[i],gnlist_com_dwgene_6v12$source[i],
                                                      gnlist_com_dwgene_6v12$p_value[i],gnlist_com_dwgene_6v12$term_id[i],
                                                      gnlist_com_dwgene_6v12$term_name[i],gnlist_com_dwgene_6v12$query_size[i]
                                                      ,gnlist_com_dwgene_6v12$intersection_size[i],gnlist_com_dwgene_6v12$intersection[i]))
  }
}

colnames(gnlist_com_6v12h_dw2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_6v12h_dw2$p_value=round(as.numeric(gnlist_com_6v12h_dw2$p_value) , digits = 5)
gnlist_com_6v12h_dw2=gnlist_com_6v12h_dw2[order(as.numeric(gnlist_com_6v12h_dw2$p_value)),]
write.csv(gnlist_com_6v12h_dw2, file = "gnlist_com_6v12h_dw2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# 6v24 dw  too lass
gnlist_com_6v24h_dw2<- data.frame()
for (i in 1:length(gnlist_com_dwgene_6v24$term_size)){
  if(gnlist_com_dwgene_6v24$term_size[i]<1000 && gnlist_com_dwgene_6v24$term_size[i]>10){
    print(i)
    print(gnlist_com_dwgene_6v24$term_size[i])
    gnlist_com_6v24h_dw2=rbind(gnlist_com_6v24h_dw2,c(gnlist_com_dwgene_6v24$term_size[i],gnlist_com_dwgene_6v24$source[i],
                                                      gnlist_com_dwgene_6v24$p_value[i],gnlist_com_dwgene_6v24$term_id[i],
                                                      gnlist_com_dwgene_6v24$term_name[i],gnlist_com_dwgene_6v24$query_size[i]
                                                      ,gnlist_com_dwgene_6v24$intersection_size[i],gnlist_com_dwgene_6v24$intersection[i]))
  }
}

colnames(gnlist_com_6v24h_dw2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_6v24h_dw2$p_value=round(as.numeric(gnlist_com_6v24h_dw2$p_value) , digits = 5)
gnlist_com_6v24h_dw2=gnlist_com_6v24h_dw2[order(as.numeric(gnlist_com_6v24h_dw2$p_value)),]
write.csv(gnlist_com_6v24h_dw2, file = "gnlist_com_6v24h_dw2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#

gnlist_com_12v24h_dw2<- data.frame()
for (i in 1:length(gnlist_com_dwgene_12v24$term_size)){
  if(gnlist_com_dwgene_12v24$term_size[i]<1000 && gnlist_com_dwgene_12v24$term_size[i]>10){
    print(i)
    print(gnlist_com_dwgene_12v24$term_size[i])
    gnlist_com_12v24h_dw2=rbind(gnlist_com_12v24h_dw2,c(gnlist_com_dwgene_12v24$term_size[i],gnlist_com_dwgene_12v24$source[i],
                                                        gnlist_com_dwgene_12v24$p_value[i],gnlist_com_dwgene_12v24$term_id[i],
                                                        gnlist_com_dwgene_12v24$term_name[i],gnlist_com_dwgene_12v24$query_size[i]
                                                        ,gnlist_com_dwgene_12v24$intersection_size[i],gnlist_com_dwgene_12v24$intersection[i]))
  }
}

colnames(gnlist_com_12v24h_dw2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_12v24h_dw2$p_value=round(as.numeric(gnlist_com_12v24h_dw2$p_value) , digits = 5)
gnlist_com_12v24h_dw2=gnlist_com_12v24h_dw2[order(as.numeric(gnlist_com_12v24h_dw2$p_value)),]
write.csv(gnlist_com_12v24h_dw2, file = "gnlist_com_12v24h_dw2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

gnlist_com_6v12v24h_dw2<- data.frame()
for (i in 1:length(gnlist_com_dwgene_6v12v24$term_size)){
  if(gnlist_com_dwgene_6v12v24$term_size[i]<1000 && gnlist_com_dwgene_6v12v24$term_size[i]>10){
    print(i)
    print(gnlist_com_dwgene_6v12v24$term_size[i])
    gnlist_com_6v12v24h_dw2=rbind(gnlist_com_6v12v24h_dw2,c(gnlist_com_dwgene_6v12v24$term_size[i],gnlist_com_dwgene_6v12v24$source[i],
                                                            gnlist_com_dwgene_6v12v24$p_value[i],gnlist_com_dwgene_6v12v24$term_id[i],
                                                            gnlist_com_dwgene_6v12v24$term_name[i],gnlist_com_dwgene_6v12v24$query_size[i]
                                                            ,gnlist_com_dwgene_6v12v24$intersection_size[i],gnlist_com_dwgene_6v12v24$intersection[i]))
  }
}

colnames(gnlist_com_6v12v24h_dw2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
gnlist_com_6v12v24h_dw2$p_value=round(as.numeric(gnlist_com_6v12v24h_dw2$p_value) , digits = 5)
gnlist_com_6v12v24h_dw2=gnlist_com_6v12v24h_dw2[order(as.numeric(gnlist_com_6v12v24h_dw2$p_value)),]
write.csv(gnlist_com_6v12v24h_dw2, file = "gnlist_com_6v12v24h_dw2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


# change to character to write csv  to save total geneset
gnlist_gnname_6h_up2$parents= as.character(gnlist_gnname_6h_up2$parents) 
gnlist_gnname_12h_up2$parents= as.character(gnlist_gnname_12h_up2$parents) 
gnlist_gnname_24h_up2$parents= as.character(gnlist_gnname_24h_up2$parents)

gnlist_com_upgene_6v12$parents= as.character(gnlist_com_6v12h_up2$parents) 
gnlist_com_upgene_6v24$parents= as.character(gnlist_com_upgene_6v24$parents) 
gnlist_com_upgene_12v24$parents= as.character(gnlist_com_upgene_12v24$parents) 
gnlist_com_upgene_6v12v24$parents= as.character(gnlist_com_upgene_6v12v24$parents) 
#dw
gnlist_gnname_6h_dw2$parents= as.character(gnlist_gnname_6h_dw2$parents) 
gnlist_gnname_12h_dw2$parents= as.character(gnlist_gnname_12h_dw2$parents) 
gnlist_gnname_24h_dw2$parents= as.character(gnlist_gnname_24h_dw2$parents)

gnlist_com_dwgene_6v12$parents= as.character(gnlist_com_dwgene_6v12$parents) 
gnlist_com_dwgene_6v24$parents= as.character(gnlist_com_dwgene_6v24$parents) 
gnlist_com_dwgene_12v24$parents= as.character(gnlist_com_dwgene_12v24$parents) 
gnlist_com_dwgene_6v12v24$parents= as.character(gnlist_com_dwgene_6v12v24$parents) 

# save table with total
write.csv(gnlist_gnname_6h_up, file = "genelist_genename_6h_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_12h_up, file = "genelist_genename_12h_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_24h_up, file = "genelist_genename_24h_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)

write.csv(gnlist_com_upgene_6v12, file = "genelist_com_6v12_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_upgene_12v24, file = "genelist_com_12v24_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_upgene_6v24, file = "genelist_com_6v24_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_upgene_6v12v24, file = "genelist_com_6v12v24_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
# save table  dw
write.csv(gnlist_gnname_6h_dw, file = "genelist_genename_6h_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_12h_dw, file = "genelist_genename_12h_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_24h_dw, file = "genelist_genename_24h_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)

write.csv(gnlist_com_dwgene_6v12, file = "genelist_com_6v12_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_dwgene_6v24, file = "genelist_com_6v24_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_dwgene_12v24, file = "genelist_com_12v24_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_dwgene_6v12v24, file = "genelist_com_6v12v24_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)


# dot plot preprocess
dotp_6h_up=data.frame(gnlist_gnname_6h_up2$term_name,gnlist_gnname_6h_up2$p_value) 
dotp_12h_up=data.frame(gnlist_gnname_12h_up2$term_name,gnlist_gnname_12h_up2$p_value) 
dotp_24h_up=data.frame(gnlist_gnname_24h_up2$term_name,gnlist_gnname_24h_up2$p_value) 

dotp_6v12_up=data.frame(gnlist_com_6v12h_up2$term_name,gnlist_com_6v12h_up2$p_value) 
dotp_6v24_up=data.frame(gnlist_com_6v24h_up2$term_name,gnlist_com_6v24h_up2$p_value) 
dotp_12v24_up=data.frame(gnlist_com_12v24h_up2$term_name,gnlist_com_12v24h_up2$p_value) 
dotp_6v12v24_up=data.frame(gnlist_com_6v12v24h_up2$term_name,gnlist_com_6v12v24h_up2$p_value) 

# dot plot preprocess  down regulate
dotp_6h_dw=data.frame(gnlist_gnname_6h_dw2$term_name,gnlist_gnname_6h_dw2$p_value) 
dotp_12h_dw=data.frame(gnlist_gnname_12h_dw2$term_name,gnlist_gnname_12h_dw2$p_value) 
dotp_24h_dw=data.frame(gnlist_gnname_24h_dw2$term_name,gnlist_gnname_24h_dw2$p_value) 

dotp_6v12_dw=data.frame(gnlist_com_6v12h_dw2$term_name,gnlist_com_6v12h_dw2$p_value) 
dotp_6v24_dw=data.frame(gnlist_com_6v24h_dw2$term_name,gnlist_com_6v24h_dw2$p_value)  # no vector
dotp_12v24_dw=data.frame(gnlist_com_12v24h_dw2$term_name,gnlist_com_12v24h_dw2$p_value) 
dotp_6v12v24_dw=data.frame(gnlist_com_6v12v24h_dw2$term_name,gnlist_com_6v12v24h_dw2$p_value) 

# order data by p value
dotp_6h_up=dotp_6h_up[order(dotp_6h_up$gnlist_gnname_6h_up2.p_value),]
dotp_12h_up=dotp_12h_up[order(dotp_12h_up$gnlist_gnname_12h_up2.p_value),]
dotp_24h_up=dotp_24h_up[order(dotp_24h_up$gnlist_gnname_24h_up2.p_value),]

dotp_6v12_up=dotp_6v12_up[order(dotp_6v12_up$gnlist_com_6v12h_up2.p_value),]
dotp_6v24_up=dotp_6v24_up[order(dotp_6v24_up$gnlist_com_6v24h_up2.p_value),]
dotp_12v24_up=dotp_12v24_up[order(dotp_12v24_up$gnlist_com_12v24h_up2.p_value),]
dotp_6v12v24_up=dotp_6v12v24_up[order(dotp_6v12v24_up$gnlist_com_6v12v24h_up2.p_value),]

# order data by p value  down
dotp_6h_dw=dotp_6h_dw[order(dotp_6h_dw$gnlist_gnname_6h_dw2.p_value),]
dotp_12h_dw=dotp_12h_dw[order(dotp_12h_dw$gnlist_gnname_12h_dw2.p_value),]
dotp_24h_dw=dotp_24h_dw[order(dotp_24h_dw$gnlist_gnname_24h_dw2.p_value),]

dotp_6v12_dw=dotp_6v12_dw[order(dotp_6v12_dw$gnlist_com_6v12h_dw2.p_value),]
dotp_6v24_dw=dotp_6v24_dw[order(dotp_6v24_dw$gnlist_com_dwgene_6v24.p_value),]
dotp_12v24_dw=dotp_12v24_dw[order(dotp_12v24_dw$gnlist_com_12v24h_dw2.p_value),]
dotp_6v12v24_dw=dotp_6v12v24_dw[order(dotp_6v12v24_dw$gnlist_com_6v12v24h_dw2.p_value),]

write.csv(dotp_6h_up, file = "genelist_pvalue_6_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12h_up, file = "genelist_pvalue_12_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_24h_up, file = "genelist_pvalue_24_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6h_dw, file = "genelist_pvalue_6_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12h_dw, file = "genelist_pvalue_12_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_24h_dw, file = "genelist_pvalue_24_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)

write.csv(dotp_6v12_up, file = "genelist_pvalue_6v12_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v24_up, file = "genelist_pvalue_6v24_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12v24_up, file = "genelist_pvalue_12v24_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v12v24_up, file = "genelist_pvalue_6v12v24_up.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v12_dw, file = "genelist_pvalue_6v12_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v24_dw, file = "genelist_pvalue_6v24_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12v24_dw, file = "genelist_pvalue_12v24_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v12v24_dw, file = "genelist_pvalue_6v12v24_dw.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)

# dot plot up  initial without order
ggplot(dotp_6h_up, aes(x = dotp_6h_up[,2], y = dotp_6h_up[,1]))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#FFAAD4", stackratio = 1.5, dotsize = 0.5)+
  labs(title= "dotplot_6h_up", x = "p value", y = NULL)

# dot plot with order
ggplot(dotp_6h_up, aes(x = dotp_6h_up[,2], y = reorder(dotp_6h_up[,1],-as.numeric(dotp_6h_up[,2])))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#FFAAD4", stackratio = 1.5, dotsize = 0.3)+
  labs(title= "dotplot_6h_up", x = "p value", y = NULL)
ggplot(dotp_12h_up, aes(x = dotp_12h_up[,2], y = reorder(dotp_12h_up[,1],-as.numeric(dotp_12h_up[,2]))))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#aafbff", stackratio = 1.5, dotsize = 0.1)+
  labs(title= "dotplot_12h_up", x = "p value", y = NULL)
ggplot(dotp_24h_up, aes(x = dotp_24h_up[,2], y = reorder(dotp_24h_up[,1],-as.numeric(dotp_24h_up[,2]))))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#ffe1aa", stackratio = 1, dotsize = 0.1)+
  labs(title= "dotplot_24h_up", x = "p value", y = NULL)

# dot plot with order down
ggplot(dotp_6h_dw, aes(x = dotp_6h_dw[,2], y = reorder(dotp_6h_dw[,1],-as.numeric(dotp_6h_dw[,2])))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#FFAAD4", stackratio = 1.5, dotsize = 0.1)+
  labs(title= "dotplot_6h_dw", x = "p value", y = NULL)
ggplot(dotp_12h_dw, aes(x = dotp_12h_dw[,2], y = reorder(dotp_12h_dw[,1],-as.numeric(dotp_12h_dw[,2])))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#aafbff", stackratio = 1.5, dotsize = 0.1)+
  labs(title= "dotplot_12h_dw", x = "p value", y = NULL)
ggplot(dotp_24h_dw, aes(x = dotp_24h_dw[,2], y = reorder(dotp_24h_dw[,1],-as.numeric(dotp_24h_dw[,2]))))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#ffe1aa", stackratio = 1, dotsize = 0.1)+
  labs(title= "dotplot_24h_dw", x = "p value", y = NULL)

# save plot  
pdf("dotplot_6h_up.pdf")
pdf("dotplot_12h_up.pdf")
pdf("dotplot_24h_up.pdf")
pdf("dotplot_6h_dw.pdf")
pdf("dotplot_12h_dw.pdf")
pdf("dotplot_24h_dw.pdf")
dev.off()

# filter lower number
dotp_6h_up2= dotp_6h_up[dotp_6h_up$gnlist_gnname_6h_up.p.value<0.01,]
dotp_12h_up2= dotp_12h_up[dotp_12h_up$gnlist_gnname_12h_up.p.value<0.01,]
dotp_24h_up2= dotp_24h_up[dotp_24h_up$gnlist_gnname_24h_up.p.value<0.01,]
dotp_6h_dw2= dotp_6h_dw[dotp_6h_dw$gnlist_gnname_6h_dw.p.value<0.01,]
dotp_12h_dw2= dotp_12h_dw[dotp_12h_dw$gnlist_gnname_12h_dw.p.value<0.01,]
dotp_24h_dw2= dotp_24h_dw[dotp_24h_dw$gnlist_gnname_24h_dw.p.value<0.01,]

#export genenames from every timepoint
mRNA_genenames= c()
mRNA_genenames= c(genenames_6h_up,genenames_12h_up,genenames_24h_up,genenames_6h_dw,genenames_12h_dw,genenames_24h_dw)
lapply(mRNA_genenames, function(x) write.table( data.frame(x), 'genenames_mRNA.csv'  , append= T, sep=',' ))

# export geneset term ID from every timepoint
gnlist_6h_up=as.data.frame(cbind(gnlist_gnname_6h_up2$term_ID, gnlist_gnname_6h_up2$term_name))
gnlist_12h_up=as.data.frame(cbind(gnlist_gnname_12h_up2$term_ID, gnlist_gnname_12h_up2$term_name))
gnlist_24h_up=as.data.frame(cbind(gnlist_gnname_24h_up2$term_ID, gnlist_gnname_24h_up2$term_name))
gnlist_6h_dw=as.data.frame(cbind(gnlist_gnname_6h_dw2$term_ID, gnlist_gnname_6h_dw2$term_name))
gnlist_12h_dw=as.data.frame(cbind(gnlist_gnname_12h_dw2$term_ID, gnlist_gnname_12h_dw2$term_name))
gnlist_24h_dw=as.data.frame(cbind(gnlist_gnname_24h_dw2$term_ID, gnlist_gnname_24h_dw2$term_name))
mRNA_gnlist= c(gnlist_6h_up,gnlist_12h_up,gnlist_24h_up,gnlist_6h_dw, gnlist_12h_dw, gnlist_24h_dw)
lapply(mRNA_gnlist, function(x) write.table( data.frame(x), 'geneset_mRNA.csv'  , append= T, sep=',' ))


# export up dw gn table

write.csv(upregulated6h, file = "up-6h-mrna.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(upregulated12h, file = "up-12h-mrna.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(upregulated24h, file = "up-24h-mrna.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(downregulated6h, file = "dw-6h-mrna.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(downregulated12h, file = "dw-12h-mrna.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(downregulated24h, file = "dw-24h-mrna.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)


# retrieve sequence of up dw ENSG
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", GRCh = 37)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>
#這個也可以 但getbm 可以拿更多  getbm序列這裏只能一次一個
seq_up_6h = getSequence(id=upregulated6h$ID, 
                        type = "ensembl_gene_id", 
                        seqType="3utr",
                        mart = mart)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>

seq_up_6h = getBM(attributes = c('ensembl_transcript_id',
                     'ensembl_gene_id',
                     '3utr'),
      filters = 'ensembl_gene_id',
      values = upregulated6h$ID,
      mart = ensembl)

seq_up_12h = getBM(attributes = c('ensembl_transcript_id',
                                 'ensembl_gene_id',
                                 '3utr'),
                  filters = 'ensembl_gene_id',
                  values = upregulated12h$ID,
                  mart = ensembl)

seq_up_24h = getBM(attributes = c('ensembl_transcript_id',
                                 'ensembl_gene_id',
                                 '3utr'),
                  filters = 'ensembl_gene_id',
                  values = upregulated24h$ID,
                  mart = ensembl)

# dw
seq_dw_6h = getBM(attributes = c('ensembl_transcript_id',
                                 'ensembl_gene_id',
                                 '3utr'),
                  filters = 'ensembl_gene_id',
                  values = downregulated6h$ID,
                  mart = ensembl)

seq_dw_12h = getBM(attributes = c('ensembl_transcript_id',
                                 'ensembl_gene_id',
                                 '3utr'),
                  filters = 'ensembl_gene_id',
                  values = downregulated12h$ID,
                  mart = ensembl)

seq_dw_24h = getBM(attributes = c('ensembl_transcript_id',
                                 'ensembl_gene_id',
                                 '3utr'),
                  filters = 'ensembl_gene_id',
                  values = downregulated24h$ID,
                  mart = ensembl)


# relocate ID
seq_up_6h <- seq_up_6h %>% relocate(ensembl_transcript_id, .after = ensembl_gene_id)
seq_up_12h <- seq_up_12h %>% relocate(ensembl_transcript_id, .after = ensembl_gene_id)
seq_up_24h <- seq_up_24h %>% relocate(ensembl_transcript_id, .after = ensembl_gene_id)

seq_dw_6h <- seq_dw_6h %>% relocate(ensembl_transcript_id, .after = ensembl_gene_id)
seq_dw_12h <- seq_dw_12h %>% relocate(ensembl_transcript_id, .after = ensembl_gene_id)
seq_dw_24h <- seq_dw_24h %>% relocate(ensembl_transcript_id, .after = ensembl_gene_id)

#export fasta
exportFASTA(seq_up_6h,file="seq_up_6h_a549_mrna.fasta")
exportFASTA(seq_up_12h,file="seq_up_12h_a549_mrna.fasta")
exportFASTA(seq_up_24h,file="seq_up_24h_a549_mrna.fasta")

exportFASTA(seq_dw_6h,file="seq_dw_6h_a549_mrna.fasta")
exportFASTA(seq_dw_12h,file="seq_dw_12h_a549_mrna.fasta")
exportFASTA(seq_dw_24h,file="seq_dw_24h_a549_mrna.fasta")
