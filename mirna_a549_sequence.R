library(BiocManager)
library(DESeq2)
library(biomaRt)

setwd("/home/yuanpeng/yuan/R sarsCoV2/mirna_sequen")

fr <- read.table(file = "A549_hg19_miRNA_only.txt", header = TRUE, sep = ",")

# change remove rowname
myrownames=fr[,1]
fr2=fr
fr2=fr2[,-c(1:6)]
rownames(fr2)=myrownames

# change col name
rsmID <- c("0_r2_notinf","0_r3_notinf","3_r2_notinf","3_r3_notinf","6_r2_notinf","6_r3_notinf","12_r2_notinf","12_r3_notinf","24_r2_notinf","24_r3_notinf","3_r2_inf","3_r3_inf","6_r2_inf","6_r3_inf","12_r2_inf","12_r3_inf","24_r2_inf","24_r3_inf")
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
condition=c("controls","controls","controls","controls","controls","controls","controls","controls","controls","controls","infected","infected","infected","infected","infected","infected","infected","infected")
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
Mnormalized3= Mnormalized3[,!colnames(Mnormalized3) %in% c("3_r2_notinf","3_r3_notinf","6_r2_notinf","6_r3_notinf","12_r2_notinf","12_r3_notinf","24_r2_notinf","24_r3_notinf")]


# clustering 
method <- "average"
distance <- "manhattan"
clustering <- hclust(dist(t(Mnormalized3),method = distance),method = method)
plot(clustering)


#remove some outlier first
Mnormalized3= Mnormalized3[,!colnames(Mnormalized3) %in% c( "3_r2_inf", "12_r2_inf")]
#cluster again
clustering <- hclust(dist(t(Mnormalized3),method = distance),method = method)
plot(clustering)

# pvclust
library(pvclust)
method <- "average"
distance <- "correlation"
clustering <- pvclust((Mnormalized3),method.hclust="average", method.dist="correlation")
plot(clustering)


# example scatterplot of normlaized data
plot(Mnormalized3[,1],Mnormalized3[,2],pch=".")

# remove outlier from raw data so can do deseq
fr3=fr2
fr3= fr3[,!colnames(fr3) %in% c("3_r2_notinf","3_r3_notinf","6_r2_notinf","6_r3_notinf","12_r2_notinf","12_r3_notinf","24_r2_notinf","24_r3_notinf", "3_r2_inf", "12_r2_inf")]
fr3=as.matrix(fr3)


#>>>>>>>>>>>> for perason corelation >>>>>>>>>>>>>
mean1=0
for(i in 1:length(fr3[,1])){
  mean1[i]=median(fr3[i,])}
fr4=fr3[mean1>quantile(mean1,0.3),]
write.csv(fr3, file = "mirna_fr3.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
#>>>>>>>>>>>>>

# 6 hrs vs control
cfr3=fr3[,c(1,2,3,4,5)]
condition=c("controls","controls","infected","infected","infected")

#Create 6hr DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(cfr3, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
norm_count = counts(dds, normalized=T)

# norm_count = as.matrix(fr3[,c(1:5)])

#0 hr sum result from deseq
sum_array_0h= c()
for (row in 1:length(norm_count[,1])){
  sum_array_0h= c(sum_array_0h, sum(c(norm_count[row,1],
                                      norm_count[row,2])))
}
sum_sum_0h= sum(sum_array_0h)
# 6hr sum result from deseq
sum_array_6h= c()
for (row in 1:length(norm_count[,1])){
  sum_array_6h= c(sum_array_6h, sum(c(norm_count[row,5],norm_count[row,4],
                                      norm_count[row,3])))
}
sum_sum_6h= sum(sum_array_6h)

# create table for chi test
# 0 vs 6hr
sum_result_6h=c()
for (row in 1:length(sum_array_0h)) {
  col_1= c(sum_array_0h[row],sum_sum_0h-sum_array_0h[row])
  col_2= c(sum_array_6h[row], sum_sum_6h-sum_array_6h[row])
  matrix= as.matrix(c(col_1,col_2))
  dim(matrix)=c(2,2)
  temp_row=c(chisq.test(matrix)$p.value,sum_array_0h[row],sum_array_6h[row])
  sum_result_6h=rbind(sum_result_6h,temp_row)
}

sum_result_6h=as.matrix(sum_result_6h)
row.names(sum_result_6h)= row.names(norm_count)
colnames(sum_result_6h)=c("pvalue","norm_counts_0h","norm_counts_6h")

# order p value from small to big  and keep p <0.05
sum_result_6h=sum_result_6h[order(sum_result_6h[,1]), ,drop=F]
sum_result_6h2=sum_result_6h[sum_result_6h[,1]<=0.05, ,drop=F]

# add log2 fold change
FC=c()
for (row in 1:length(sum_result_6h[,1])) {
  FC[row]=as.numeric(log2((sum_result_6h[row,3]/sum(sum_result_6h[,3]))/(sum_result_6h[row,2]/sum(sum_result_6h[,2]))))
  
}
sum_result_6h= cbind(sum_result_6h,as.data.frame(FC))
colnames(sum_result_6h)=c("pvalue","norm_counts_0h","norm_counts_6h","log2FoldChange")

FC2=c()
for (row in 1:length(sum_result_6h2[,1])) {
  FC2[row]=as.numeric(log2((sum_result_6h2[row,3]/sum(sum_result_6h2[,3]))/(sum_result_6h2[row,2]/sum(sum_result_6h2[,2]))))
  
}
sum_result_6h2= cbind(sum_result_6h2,as.data.frame(FC2))
colnames(sum_result_6h2)=c("pvalue","norm_counts_0h","norm_counts_6h","log2FoldChange")


# select up down regulate
upregulated6h=c()
downregulated6h=c()
for (row in 1:length(sum_result_6h2[,1])) {
  if (as.numeric(sum_result_6h2[row,4]) > 0.05) {
    upregulated6h= rbind(upregulated6h,sum_result_6h2[row,,drop=F])
  }  
  if (as.numeric(sum_result_6h[row,4]) < -0.05) {
    downregulated6h=rbind(downregulated6h,sum_result_6h2[row,,drop=F])
  }
}

#volcano plot total
temp_res= sum_result_6h
temp_res$diffexpressed= "No"
temp_res$diffexpressed[temp_res$log2FoldChange>0.05 & temp_res$pvalue<0.05]="Up"
temp_res$diffexpressed[temp_res$log2FoldChange< -0.05 & temp_res$pvalue<0.05]="Down"
temp_res_vol= ggplot(temp_res, aes(x=temp_res$log2FoldChange, y=-log10(temp_res$pvalue)))+ geom_point() + theme_minimal()
temp_res_vol2= temp_res_vol+ geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_res_vol3 = temp_res_vol2 + scale_color_manual(values=mycolors)

ggplot(temp_res, aes(x=temp_res_vol3$data$log2FoldChange, y=-log10(temp_res_vol3$data$pvalue), col=temp_res_vol3$data$diffexpressed))+ 
  geom_point() + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")+ 
  theme_minimal()+labs(title= "Volcano_plot_6h", x = "log2FoldChange", y = "-log10(pvalue)", colour="Diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))

table(temp_res$diffexpressed)

#>>>>>>>>>>>>>>>>>>>>>>>>
#volcano plot
temp_6_updw=rbind(upregulated6h,downregulated6h)
temp_6_updw_vol= ggplot(temp_6_updw, aes(x=temp_6_updw$log2FoldChange, y=-log10(temp_6_updw$pvalue)))+ geom_point() + theme_minimal()
temp_6_updw_vol2 = temp_6_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

temp_6_updw$diffexpressed= "No"
temp_6_updw$diffexpressed[temp_6_updw$log2FoldChange>0.05 & temp_6_updw$pvalue<0.05]="Up"
temp_6_updw$diffexpressed[temp_6_updw$log2FoldChange<0.05 & temp_6_updw$pvalue<0.05]="Down"

temp_6_updw_vol= ggplot(temp_6_updw, aes(x=temp_6_updw$log2FoldChange, y=-log10(temp_6_updw$pvalue), col=temp_6_updw$diffexpressed))+ geom_point() + theme_minimal()
temp_6_updw_vol2 = temp_6_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_6_updw_vol3=temp_6_updw_vol2+ scale_color_manual(values=mycolors)

ggplot(temp_6_updw, aes(x=temp_6_updw_vol3$data$log2FoldChange, y=-log10(temp_6_updw_vol3$data$pvalue), col=temp_6_updw_vol3$data$diffexpressed))+ 
  geom_point() + theme_minimal()+labs(title= "volcano_plot_6h", x = "log2FoldChange", y = "-log10(pvalue)", colour="diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+ 
  scale_x_continuous(limits=c(-5,5), breaks=seq(-5,5,2) )
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



# 12 hrs vs control
cfr3=fr3[,c(1,2,6)]
condition=c("controls","controls","infected")

#Create 12hr DESeq2Dataset object 
dds <- DESeqDataSetFromMatrix(cfr3, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
norm_count = counts(dds, normalized=T)

# norm_count = as.matrix(fr3[,c(1,2,6)])

#0 hr sum result from deseq
sum_array_0h= c()
for (row in 1:length(norm_count[,1])){
  sum_array_0h= c(sum_array_0h, sum(c(norm_count[row,1],
                                      norm_count[row,2])))
}
sum_sum_0h= sum(sum_array_0h)


# 12hr sum result from deseq
sum_array_12h= c()
for (row in 1:length(norm_count[,1])){
  sum_array_12h= c(sum_array_12h, sum(c(norm_count[row,3])))
}
sum_sum_12h= sum(sum_array_12h)

# 0 vs 12hr
sum_result_12h=c()
for (row in 1:length(sum_array_0h)) {
  col_1= c(sum_array_0h[row],sum_sum_0h-sum_array_0h[row])
  col_2= c(sum_array_12h[row], sum_sum_12h-sum_array_12h[row])
  matrix= as.matrix(c(col_1,col_2))
  dim(matrix)=c(2,2)
  temp_row=c(chisq.test(matrix)$p.value,sum_array_0h[row],sum_array_12h[row])
  sum_result_12h=rbind(sum_result_12h,temp_row)
}

sum_result_12h=as.matrix(sum_result_12h)
row.names(sum_result_12h)= row.names(norm_count)
colnames(sum_result_12h)=c("pvalue","norm_counts_0h","norm_counts_12h")

# order p value from small to big  and keep <0.05
sum_result_12h=sum_result_12h[order(sum_result_12h[,1]), ,drop=F]
sum_result_12h2=sum_result_12h[sum_result_12h[,1]<0.05, ,drop=F]

# add log2 fold change
FC=c()
for (row in 1:length(sum_result_12h[,1])) {
  FC[row]=as.numeric(log2((sum_result_12h[row,3]/sum(sum_result_12h[,3]))/(sum_result_12h[row,2]/sum(sum_result_12h[,2]))))
  
}
sum_result_12h= cbind(sum_result_12h,as.data.frame(FC))
colnames(sum_result_12h)=c("pvalue","norm_counts_0h","norm_counts_12h","log2FoldChange")

FC2=c()
for (row in 1:length(sum_result_12h2[,1])) {
  FC2[row]=as.numeric(log2((sum_result_12h[row,3]/sum(sum_result_12h[,3]))/(sum_result_12h[row,2]/sum(sum_result_12h[,2]))))
  
}
sum_result_12h2= cbind(sum_result_12h2,as.data.frame(FC2))
colnames(sum_result_12h2)=c("pvalue","norm_counts_0h","norm_counts_12h","log2FoldChange")

# select up down regulate
upregulated12h=c()
downregulated12h=c()
for (row in 1:length(sum_result_12h2[,1])) {
  if (as.numeric(sum_result_12h2[row,4]) > 0.05) {
    upregulated12h= rbind(upregulated12h,sum_result_12h2[row,,drop=F])
  }  
  if (as.numeric(sum_result_12h[row,4]) < -0.05) {
    downregulated12h=rbind(downregulated12h,sum_result_12h2[row,,drop=F])
  }
}

#volcano plot total
temp_res= sum_result_12h
temp_res$diffexpressed= "No"
temp_res$diffexpressed[temp_res$log2FoldChange>0.05 & temp_res$pvalue<0.05]="Up"
temp_res$diffexpressed[temp_res$log2FoldChange< -0.05 & temp_res$pvalue<0.05]="Down"

temp_res_vol= ggplot(temp_res, aes(x=temp_res$log2FoldChange, y=-log10(temp_res$pvalue)))+ geom_point() + theme_minimal()
temp_res_vol2= temp_res_vol+ geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_res_vol3 = temp_res_vol2 + scale_color_manual(values=mycolors)

ggplot(temp_res, aes(x=temp_res_vol3$data$log2FoldChange, y=-log10(temp_res_vol3$data$pvalue), col=temp_res_vol3$data$diffexpressed))+ 
  geom_point() + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")+ 
  theme_minimal()+labs(title= "Volcano_plot_12h", x = "log2FoldChange", y = "-log10(pvalue)", colour="Diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))

table(temp_res$diffexpressed)

#>>>>>>>>>>>>>>>>>>>>>
#volcano plot
temp_12_updw=rbind(upregulated12h,downregulated12h)
temp_12_updw_vol= ggplot(temp_12_updw, aes(x=temp_12_updw$log2FoldChange, y=-log10(temp_12_updw$pvalue)))+ geom_point() + theme_minimal()
temp_12_updw_vol2 = temp_12_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

temp_12_updw$diffexpressed= "NO"
temp_12_updw$diffexpressed[temp_12_updw$log2FoldChange>0.05 & temp_12_updw$pvalue<0.05]="up"
temp_12_updw$diffexpressed[temp_12_updw$log2FoldChange<0.05 & temp_12_updw$pvalue<0.05]="down"

temp_12_updw_vol= ggplot(temp_12_updw, aes(x=temp_12_updw$log2FoldChange, y=-log10(temp_12_updw$pvalue), col=temp_12_updw$diffexpressed))+ geom_point() + theme_minimal()
temp_12_updw_vol2 = temp_12_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_12_updw_vol3=temp_12_updw_vol2+ scale_color_manual(values=mycolors)

ggplot(temp_12_updw, aes(x=temp_12_updw_vol3$data$log2FoldChange, y=-log10(temp_12_updw_vol3$data$pvalue), col=temp_12_updw_vol3$data$diffexpressed))+ 
  geom_point() + theme_minimal()+labs(title= "volcano_plot_12h", x = "log2FoldChange", y = "-log10(pvalue)", colour="diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-5,5), breaks=seq(-5,5,2) )

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# 24 hrs vs control
cfr3=fr3[,c(1,2,7,8)]
condition=c("controls","controls","infected","infected")

#Create 24hr DESeq2Dataset object  
dds <- DESeqDataSetFromMatrix(cfr3, DataFrame(condition), ~ condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
norm_count = counts(dds, normalized=T)

# norm_count = as.matrix(fr3[,c(1,2,7,8)])

#0 hr sum result from deseq
sum_array_0h= c()
for (row in 1:length(norm_count[,1])){
  sum_array_0h= c(sum_array_0h, sum(c(norm_count[row,1],
                                      norm_count[row,2])))
}
sum_sum_0h= sum(sum_array_0h)

# 24hr sum result from deseq
sum_array_24h= c()
for (row in 1:length(norm_count[,1])){
  sum_array_24h= c(sum_array_24h, sum(c(norm_count[row,3],
                                        norm_count[row,4])))
}
sum_sum_24h= sum(sum_array_24h)

# 0 vs 24hr
sum_result_24h=c()
for (row in 1:length(sum_array_0h)) {
  col_1= c(sum_array_0h[row],sum_sum_0h-sum_array_0h[row])
  col_2= c(sum_array_24h[row], sum_sum_24h-sum_array_24h[row])
  matrix= as.matrix(c(col_1,col_2))
  dim(matrix)=c(2,2)
  temp_row=c(chisq.test(matrix)$p.value,sum_array_0h[row],sum_array_24h[row])
  sum_result_24h=rbind(sum_result_24h,temp_row)
}

sum_result_24h=as.matrix(sum_result_24h)
row.names(sum_result_24h)= row.names(norm_count)
colnames(sum_result_24h)=c("pvalue","norm_counts_0h","norm_counts_24h")

# order p value from small to big  and remove <0.05
sum_result_24h=sum_result_24h[order(sum_result_24h[,1]), ,drop=F]
sum_result_24h2=sum_result_24h[sum_result_24h[,1]<0.05, ,drop=F]

# add log2 fold change
FC=c()
for (row in 1:length(sum_result_24h[,1])) {
  FC[row]=as.numeric(log2((sum_result_24h[row,3]/sum(sum_result_24h[,3]))/(sum_result_24h[row,2]/sum(sum_result_24h[,2]))))
  
}
sum_result_24h= cbind(sum_result_24h,as.data.frame(FC))
colnames(sum_result_24h)=c("pvalue","norm_counts_0h","norm_counts_24h","log2FoldChange")

FC2=c()
for (row in 1:length(sum_result_24h2[,1])) {
  FC2[row]=as.numeric(log2((sum_result_24h2[row,3]/sum(sum_result_24h2[,3]))/(sum_result_24h2[row,2]/sum(sum_result_24h2[,2]))))
  
}
sum_result_24h2= cbind(sum_result_24h2,as.data.frame(FC2))
colnames(sum_result_24h2)=c("pvalue","norm_counts_0h","norm_counts_12h","log2FoldChange")

# select up down regulate
upregulated24h=c()
downregulated24h=c()
for (row in 1:length(sum_result_24h2[,1])) {
  if (as.numeric(sum_result_24h2[row,4]) > 0.05) {
    upregulated24h= rbind(upregulated24h,sum_result_24h2[row,,drop=F])
  }  
  if (as.numeric(sum_result_24h2[row,4]) < -0.05) {
    downregulated24h=rbind(downregulated24h,sum_result_24h2[row,,drop=F])
  }
}

#volcano plot total
temp_res= sum_result_24h
temp_res$diffexpressed= "No"
temp_res$diffexpressed[temp_res$log2FoldChange>0.05 & temp_res$pvalue<0.05]="Up"
temp_res$diffexpressed[temp_res$log2FoldChange< -0.05 & temp_res$pvalue<0.05]="Down"

temp_res_vol= ggplot(temp_res, aes(x=temp_res$log2FoldChange, y=-log10(temp_res$pvalue)))+ geom_point() + theme_minimal()
temp_res_vol2= temp_res_vol+ geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_res_vol3 = temp_res_vol2 + scale_color_manual(values=mycolors)

ggplot(temp_res, aes(x=temp_res_vol3$data$log2FoldChange, y=-log10(temp_res_vol3$data$pvalue), col=temp_res_vol3$data$diffexpressed))+ 
  geom_point() + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")+ 
  theme_minimal()+labs(title= "Volcano_plot_24h", x = "log2FoldChange", y = "-log10(pvalue)", colour="Diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-4,4), breaks=seq(-4,4,2))

table(temp_res$diffexpressed)

#>>>>>>>>>>>>>>>>>>>>>>>>>>
#volcano plot
temp_24_updw=rbind(upregulated24h,downregulated24h)
temp_24_updw_vol= ggplot(temp_24_updw, aes(x=temp_24_updw$log2FoldChange, y=-log10(temp_24_updw$pvalue)))+ geom_point() + theme_minimal()
temp_24_updw_vol2 = temp_24_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

temp_24_updw$diffexpressed= "NO"
temp_24_updw$diffexpressed[temp_24_updw$log2FoldChange>0.05 & temp_24_updw$pvalue<0.05]="up"
temp_24_updw$diffexpressed[temp_24_updw$log2FoldChange<0.05 & temp_24_updw$pvalue<0.05]="down"

temp_24_updw_vol= ggplot(temp_24_updw, aes(x=temp_24_updw$log2FoldChange, y=-log10(temp_24_updw$pvalue), col=temp_24_updw$diffexpressed))+ geom_point() + theme_minimal()
temp_24_updw_vol2 = temp_24_updw_vol + geom_vline(xintercept=c(-0.05, 0.05), col="red") +  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
temp_24_updw_vol3=temp_24_updw_vol2+ scale_color_manual(values=mycolors)

ggplot(temp_24_updw, aes(x=temp_24_updw_vol3$data$log2FoldChange, y=-log10(temp_24_updw_vol3$data$pvalue), col=temp_24_updw_vol3$data$diffexpressed))+ 
  geom_point() + theme_minimal()+labs(title= "volcano_plot_24h", x = "log2FoldChange", y = "-log10(pvalue)", colour="diffexpressed")+theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits=c(-5,5), breaks=seq(-5,5,2))


#>>>>>>>>>>>>>>>>>>

# convert ENST to ENSG
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

miRNA_ensg_6h_up=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_transcript_id_version', 
                                                'ensembl_gene_id','external_gene_name', "mirbase_id"),
                                 filters = 'ensembl_transcript_id', 
                                 values = rownames(upregulated6h),
                                 mart = mart)

miRNA_ensg_12h_up=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_transcript_id_version', 
                                                 'ensembl_gene_id','external_gene_name', "mirbase_id"),
                                  filters = 'ensembl_transcript_id', 
                                  values = rownames(upregulated12h),
                                  mart = mart)

miRNA_ensg_24h_up=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_transcript_id_version', 
                                                 'ensembl_gene_id','external_gene_name', "mirbase_id"),
                                  filters = 'ensembl_transcript_id', 
                                  values = rownames(upregulated24h),
                                  mart = mart)

miRNA_ensg_6h_dw=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_transcript_id_version', 
                                                'ensembl_gene_id','external_gene_name', "mirbase_id"),
                                 filters = 'ensembl_transcript_id', 
                                 values = rownames(downregulated6h),
                                 mart = mart)

miRNA_ensg_12h_dw=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_transcript_id_version', 
                                                 'ensembl_gene_id','external_gene_name', "mirbase_id"),
                                  filters = 'ensembl_transcript_id', 
                                  values = rownames(downregulated12h),
                                  mart = mart)

miRNA_ensg_24h_dw=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_transcript_id_version', 
                                                 'ensembl_gene_id','external_gene_name', "mirbase_id"),
                                  filters = 'ensembl_transcript_id', 
                                  values = rownames(downregulated24h),
                                  mart = mart)

# >>>>>>>>>>>>>>>>>>>拿 mir id >>>>>>>>>>>>>> 0918
mirbase_gnlist_6h_up=as.data.frame(cbind(miRNA_ensg_6h_up$ensembl_transcript_id, miRNA_ensg_6h_up$ensembl_gene_id,
                                         miRNA_ensg_6h_up$external_gene_name,miRNA_ensg_6h_up$mirbase_id))
mirbase_gnlist_12h_up=as.data.frame(cbind(miRNA_ensg_12h_up$ensembl_transcript_id, miRNA_ensg_12h_up$ensembl_gene_id,
                                          miRNA_ensg_12h_up$external_gene_name,miRNA_ensg_12h_up$mirbase_id))
mirbase_gnlist_24h_up=as.data.frame(cbind(miRNA_ensg_24h_up$ensembl_transcript_id, miRNA_ensg_24h_up$ensembl_gene_id,
                                          miRNA_ensg_24h_up$external_gene_name,miRNA_ensg_24h_up$mirbase_id))
mirbase_gnlist_6h_dw=as.data.frame(cbind(miRNA_ensg_6h_dw$ensembl_transcript_id, miRNA_ensg_6h_dw$ensembl_gene_id,
                                         miRNA_ensg_6h_dw$external_gene_name,miRNA_ensg_6h_dw$mirbase_id))
mirbase_gnlist_12h_dw=as.data.frame(cbind(miRNA_ensg_12h_dw$ensembl_transcript_id, miRNA_ensg_12h_dw$ensembl_gene_id,
                                          miRNA_ensg_12h_dw$external_gene_name,miRNA_ensg_12h_dw$mirbase_id))
mirbase_gnlist_24h_dw=as.data.frame(cbind(miRNA_ensg_24h_dw$ensembl_transcript_id, miRNA_ensg_24h_dw$ensembl_gene_id,
                                          miRNA_ensg_24h_dw$external_gene_name,miRNA_ensg_24h_dw$mirbase_id))
miRNA_mirbase_gnlist= c(mirbase_gnlist_6h_up,mirbase_gnlist_12h_up,mirbase_gnlist_24h_up,mirbase_gnlist_6h_dw, mirbase_gnlist_12h_dw, mirbase_gnlist_24h_dw)
lapply(miRNA_mirbase_gnlist, function(x) write.table( data.frame(x), 'mirbase_gnlist_miRNA2.csv'  , append= T, sep=',' ))

# >>>>>>>>>>>>>>>>>拿 mir id >>>>>>>>>>>>>>

# 0920 2023 redo >>>>>>>>>>>>>>>>>>>>>>>
r_up6h=row.names(upregulated6h)
r_up12h=row.names(upregulated12h)
r_up24h=row.names(upregulated24h)

upgene= c(r_up6h,r_up12h,r_up24h)
upgene2 <- unique(upgene)

r_dw6h=row.names(downregulated6h)
r_dw12h=row.names(downregulated12h)
r_dw24h=row.names(downregulated24h)

dwgene= c(r_dw6h,r_dw12h,r_dw24h)
dwgene2 <- unique(dwgene)

upgn_fr=fr3[row.names(fr3) %in% upgene2, ]  
dwgn_fr=fr3[row.names(fr3) %in% dwgene2, ]  

write.csv(upgn_fr, file = "upgn_fr_mi.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
write.csv(dwgn_fr, file = "dwgn_fr_mi.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)

#get gn id
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
miRNA_up=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id','external_gene_name', "mirbase_id"),
                                filters = 'ensembl_transcript_id', 
                                values = upgene2,
                                mart = ensembl)
miRNA_dw=biomaRt::getBM(attributes = c('ensembl_transcript_id', 
                                       'ensembl_gene_id','external_gene_name', "mirbase_id"),
                        filters = 'ensembl_transcript_id', 
                        values = dwgene2,
                        mart = mart)

write.csv(miRNA_up, file = "mirbase_mi_up.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
write.csv(miRNA_dw, file = "mirbase_mi_dw.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)

#正常應該可以用 但部之爲啥現在不能用
# 發現是因為r 開起來太久 temp被山了
#重跑這個  就可以了
tempdir()
dir.create(tempdir())

# 0920 2023 redo >>>>>>>>>>>>>>>>>>>>>>>



# after convert to EMSG then covert to gnname
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", GRCh = 37)

# make object to store list after it convert to gene name
genenames_6h_up = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = miRNA_ensg_6h_up$ensembl_gene_id, mart = ensembl)
genenames_12h_up = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = miRNA_ensg_12h_up$ensembl_gene_id, mart = ensembl)
genenames_24h_up = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = miRNA_ensg_24h_up$ensembl_gene_id, mart = ensembl)

genenames_6h_dw = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = miRNA_ensg_6h_dw$ensembl_gene_id, mart = ensembl)
genenames_12h_dw = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = miRNA_ensg_12h_dw$ensembl_gene_id, mart = ensembl)
genenames_24h_dw = getBM(attributes = c("hgnc_symbol"),filters = "ensembl_gene_id",values = miRNA_ensg_24h_dw$ensembl_gene_id, mart = ensembl)




# venn diagram  make list 
venn_uplist= list ("inf_6hr_up"=genenames_6h_up$hgnc_symbol,"inf_12hr_up"=genenames_12h_up$hgnc_symbol,"inf_24hr_up"=genenames_24h_up$hgnc_symbol)
venn.diagram(venn_uplist, filename = 'venn_up_ensg.png', imagetype = 'png', fill = c('red', 'blue','green'), alpha = 0.50, cat.col = c('red', 'blue', 'green'), 
             col = c('red', 'blue', 'green'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')

venn_dwlist= list ("inf_6hr_dw"=genenames_6h_dw$hgnc_symbol,"inf_12hr_dw"=genenames_12h_dw$hgnc_symbol,"inf_24hr_dw"=genenames_24h_dw$hgnc_symbol)
venn.diagram(venn_dwlist, filename = 'venn_dw_ensg2.png', imagetype = 'png', fill = c('red', 'blue','green'), alpha = 0.50, cat.col = c('red', 'blue', 'green'), 
             col = c('red', 'blue', 'green'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')

#upset plot -upregulate 
pdf("upsetplot_up_mi_ensg2.pdf")
upset(fromList(venn_uplist), mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(T,F))
dev.off()

# upset plot -dwregulate 
pdf("upsetplot_dw_mi_ensg.pdf")
upset(fromList(venn_dwlist), mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(TRUE,F))
dev.off()

#extract common gene
library(ComplexHeatmap)
comdata_up= make_comb_mat(venn_uplist)
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
# genelist from up  GO:BP search from Biological Process 
gnlist_gnname_6h_up=gprofiler2::gost(genenames_6h_up[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_12h_up=gprofiler2::gost(genenames_12h_up[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_24h_up=gprofiler2::gost(genenames_24h_up[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
# only receive $result from gprofiler
gnlist_gnname_6h_up= gnlist_gnname_6h_up$result
gnlist_gnname_12h_up= gnlist_gnname_12h_up$result
gnlist_gnname_24h_up= gnlist_gnname_24h_up$result

# common up
gnlist_com_upgene_6v12= gprofiler2::gost(com_upgene_6v12, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_com_upgene_6v24= gprofiler2::gost(com_upgene_6v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)  # not significant
gnlist_com_upgene_12v24= gprofiler2::gost(com_upgene_12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)    # not significant
gnlist_com_upgene_6v12v24= gprofiler2::gost(com_upgene_6v12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)

gnlist_com_upgene_6v12= gnlist_com_upgene_6v12$result
gnlist_com_upgene_6v24= gnlist_com_upgene_6v24$result
gnlist_com_upgene_12v24= gnlist_com_upgene_12v24$result
gnlist_com_upgene_6v12v24= gnlist_com_upgene_6v12v24$result


# genelist from dw
gnlist_gnname_6h_dw=gprofiler2::gost(genenames_6h_dw[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_12h_dw=gprofiler2::gost(genenames_12h_dw[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_24h_dw=gprofiler2::gost(genenames_24h_dw[,1], correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_gnname_6h_dw= gnlist_gnname_6h_dw$result
gnlist_gnname_12h_dw= gnlist_gnname_12h_dw$result
gnlist_gnname_24h_dw= gnlist_gnname_24h_dw$result

# common dw
gnlist_com_dwgene_6v12= gprofiler2::gost(com_dwgene_6v12, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_com_dwgene_6v24= gprofiler2::gost(com_dwgene_6v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
gnlist_com_dwgene_12v24= gprofiler2::gost(com_dwgene_12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)  # not significant
gnlist_com_dwgene_6v12v24= gprofiler2::gost(com_dwgene_6v12v24, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)

gnlist_com_dwgene_6v12= gnlist_com_dwgene_6v12$result
gnlist_com_dwgene_6v24= gnlist_com_dwgene_6v24$result
gnlist_com_dwgene_12v24= gnlist_com_dwgene_12v24$result
gnlist_com_dwgene_6v12v24= gnlist_com_dwgene_6v12v24$result


# select only 10-1000 of term size
gnlist_gnname_6h_up2<- data.frame()
for (i in 1:length(gnlist_gnname_6h_up$term_size)){
  if(gnlist_gnname_6h_up$term_size[i]<1000 && gnlist_gnname_6h_up$term_size[i]>10){
    print(i)
    print(gnlist_gnname_6h_up$term_size[i])
    gnlist_gnname_6h_up2=rbind(gnlist_gnname_6h_up2,c(gnlist_gnname_6h_up$term_size[i],gnlist_gnname_6h_up$source[i],
                                                      gnlist_gnname_6h_up$p_value[i],gnlist_gnname_6h_up$term_id[i],
                                                      gnlist_gnname_6h_up$term_name[i],gnlist_gnname_6h_up$query_size[i]
                                                      ,gnlist_gnname_6h_up$intersection_size[i],gnlist_gnname_6h_up$intersection[i]))
  }
}

colnames(gnlist_gnname_6h_up2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
# order by p value
gnlist_gnname_6h_up2=gnlist_gnname_6h_up2[order(as.numeric(gnlist_gnname_6h_up2$p_value)),]
# change num size of p value to become only 5 number > 0.00001
gnlist_gnname_6h_up2$p_value=round(as.numeric(gnlist_gnname_6h_up2$p_value) , digits = 5)
write.csv(gnlist_gnname_6h_up2, file = "gnlist_gnname_6h_up2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

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
write.csv(gnlist_gnname_12h_up2, file = "gnlist_gnname_12h_up2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

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
write.csv(gnlist_gnname_24h_up2, file = "gnlist_gnname_24h_up2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#6v12
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
write.csv(gnlist_com_6v12h_up2, file = "gnlist_com_6v12h_up2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#6v24
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
gnlist_com_6v24h_up2=gnlist_com_6v24h_up2[order(as.numeric(gnlist_com_6v24h_up2$p_value)),]
gnlist_com_6v24h_up2$p_value=round(as.numeric(gnlist_com_6v24h_up2$p_value) , digits = 5)
write.csv(gnlist_com_6v24h_up2, file = "gnlist_com_6v24h_up2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#12v24
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
gnlist_com_12v24h_up2=gnlist_com_12v24h_up2[order(as.numeric(gnlist_com_12v24h_up2$p_value)),]
gnlist_com_12v24h_up2$p_value=round(as.numeric(gnlist_com_12v24h_up2$p_value) , digits = 5)
write.csv(gnlist_com_12v24h_up2, file = "gnlist_com_12v24h_up2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#6v12v24
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
gnlist_com_6v12v24h_up2=gnlist_com_6v12v24h_up2[order(as.numeric(gnlist_com_6v12v24h_up2$p_value)),]
gnlist_com_6v12v24h_up2$p_value=round(as.numeric(gnlist_com_6v12v24h_up2$p_value) , digits = 5)
write.csv(gnlist_com_6v12v24h_up2, file = "gnlist_com_6v12v24h_up2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#6 dw
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
gnlist_gnname_6h_dw2=gnlist_gnname_6h_dw2[order(as.numeric(gnlist_gnname_6h_dw2$p_value)),]
gnlist_gnname_6h_dw2$p_value=round(as.numeric(gnlist_gnname_6h_dw2$p_value) , digits = 5)
write.csv(gnlist_gnname_6h_dw2, file = "gnlist_gnname_6h_dw2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

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
gnlist_gnname_12h_dw2=gnlist_gnname_12h_dw2[order(as.numeric(gnlist_gnname_12h_dw2$p_value)),]
gnlist_gnname_12h_dw2$p_value=round(as.numeric(gnlist_gnname_12h_dw2$p_value) , digits = 5)
write.csv(gnlist_gnname_12h_dw2, file = "gnlist_gnname_12h_dw2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

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
gnlist_gnname_24h_dw2=gnlist_gnname_24h_dw2[order(as.numeric(gnlist_gnname_24h_dw2$p_value)),]
gnlist_gnname_24h_dw2$p_value=round(as.numeric(gnlist_gnname_24h_dw2$p_value) , digits = 5)
write.csv(gnlist_gnname_24h_dw2, file = "gnlist_gnname_24h_dw2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#6v12 dw
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
gnlist_com_6v12h_dw2=gnlist_com_6v12h_dw2[order(as.numeric(gnlist_com_6v12h_dw2$p_value)),]
gnlist_com_6v12h_dw2$p_value=round(as.numeric(gnlist_com_6v12h_dw2$p_value) , digits = 5)
write.csv(gnlist_com_6v12h_dw2, file = "gnlist_com_6v12h_dw2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# 6v24 dw  
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
gnlist_com_6v24h_dw2=gnlist_com_6v24h_dw2[order(as.numeric(gnlist_com_6v24h_dw2$p_value)),]
gnlist_com_6v24h_dw2$p_value=round(as.numeric(gnlist_com_6v24h_dw2$p_value) , digits = 5)
write.csv(gnlist_com_6v24h_dw2, file = "gnlist_com_6v24h_dw2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


#6v12 dw  no
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
gnlist_com_12v24h_dw2=gnlist_com_12v24h_dw2[order(as.numeric(gnlist_com_12v24h_dw2$p_value)),]
gnlist_com_12v24h_dw2$p_value=round(as.numeric(gnlist_com_12v24h_dw2$p_value) , digits = 5)
write.csv(gnlist_com_12v24h_dw2, file = "gnlist_com_12v24h_dw2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


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
gnlist_com_6v12v24h_dw2=gnlist_com_6v12v24h_dw2[order(as.numeric(gnlist_com_6v12v24h_dw2$p_value)),]
gnlist_com_6v12v24h_dw2$p_value=round(as.numeric(gnlist_com_6v12v24h_dw2$p_value) , digits = 5)
write.csv(gnlist_com_6v12v24h_dw2, file = "gnlist_com_6v12v24h_dw2_ensg.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# transfer list to dataframe character or cannot write csv
gnlist_gnname_6h_up=data.frame(lapply(gnlist_gnname_6h_up, as.character))
gnlist_gnname_12h_up=data.frame(lapply(gnlist_gnname_12h_up, as.character))
gnlist_gnname_24h_up=data.frame(lapply(gnlist_gnname_24h_up, as.character))
gnlist_gnname_6h_dw=data.frame(lapply(gnlist_gnname_6h_dw, as.character))
gnlist_gnname_12h_dw=data.frame(lapply(gnlist_gnname_12h_dw, as.character))
gnlist_gnname_24h_dw=data.frame(lapply(gnlist_gnname_24h_dw, as.character))

gnlist_com_upgene_6v12$parents= as.character(gnlist_com_upgene_6v12$parents) 
gnlist_com_upgene_6v24$parents= as.character(gnlist_com_upgene_6v24$parents) 
gnlist_com_upgene_12v24$parents= as.character(gnlist_com_upgene_12v24$parents) 
gnlist_com_upgene_6v12v24$parents= as.character(gnlist_com_upgene_6v12v24$parents) 
#dw
gnlist_com_dwgene_6v12$parents= as.character(gnlist_com_dwgene_6v12$parents) 
gnlist_com_dwgene_6v24$parents= as.character(gnlist_com_dwgene_6v24$parents) 
gnlist_com_dwgene_12v24$parents= as.character(gnlist_com_dwgene_12v24$parents) 
gnlist_com_dwgene_6v12v24$parents= as.character(gnlist_com_dwgene_6v12v24$parents)

# save table
write.table(gnlist_gnname_6h_up, file = "genelist_genename_6h_up_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_12h_up, file = "genelist_genename_12h_up_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_24h_up, file = "genelist_genename_24h_up_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
# save table  dw
write.csv(gnlist_gnname_6h_dw, file = "genelist_genename_6h_dw_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_12h_dw, file = "genelist_genename_12h_dw_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_gnname_24h_dw, file = "genelist_genename_24h_dw_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
# save table com up
write.csv(gnlist_com_upgene_6v12, file = "genelist_com_6v12_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_upgene_12v24, file = "genelist_com_12v24_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_upgene_6v24, file = "genelist_com_6v24_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_upgene_6v12v24, file = "genelist_com_6v12v24_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
# save table com dw
write.csv(gnlist_com_dwgene_6v12, file = "genelist_com_6v12_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_dwgene_6v24, file = "genelist_com_6v24_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_dwgene_12v24, file = "genelist_com_12v24_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(gnlist_com_dwgene_6v12v24, file = "genelist_com_6v12v24_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)

# dot plot preprocess
dotp_6h_up=data.frame(gnlist_gnname_6h_up2$term_name,gnlist_gnname_6h_up2$p_value) 
dotp_12h_up=data.frame(gnlist_gnname_12h_up2$term_name,gnlist_gnname_12h_up2$p_value) 
dotp_24h_up=data.frame(gnlist_gnname_24h_up2$term_name,gnlist_gnname_24h_up2$p_value) 
# dot plot preprocess  down regulate
dotp_6h_dw=data.frame(gnlist_gnname_6h_dw2$term_name,gnlist_gnname_6h_dw2$p_value) 
dotp_12h_dw=data.frame(gnlist_gnname_12h_dw2$term_name,gnlist_gnname_12h_dw2$p_value) 
dotp_24h_dw=data.frame(gnlist_gnname_24h_dw2$term_name,gnlist_gnname_24h_dw2$p_value) 
# dot com up 
dotp_6v12_up=data.frame(gnlist_com_6v12h_up2$term_name,gnlist_com_6v12h_up2$p_value) 
dotp_6v24_up=data.frame(gnlist_com_6v24h_up2$term_name,gnlist_com_6v24h_up2$p_value) 
dotp_12v24_up=data.frame(gnlist_com_12v24h_up2$term_name,gnlist_com_12v24h_up2$p_value) 
dotp_6v12v24_up=data.frame(gnlist_com_6v12v24h_up2$term_name,gnlist_com_6v12v24h_up2$p_value) 
# dot com dw 
dotp_6v12_dw=data.frame(gnlist_com_6v12h_dw2$term_name,gnlist_com_6v12h_dw2$p_value) 
dotp_6v24_dw=data.frame(gnlist_com_6v24h_dw2$term_name,gnlist_com_6v24h_dw2$p_value) 
dotp_12v24_dw=data.frame(gnlist_com_12v24h_dw2$term_name,gnlist_com_12v24h_dw2$p_value) 
dotp_6v12v24_dw=data.frame(gnlist_com_6v12v24h_dw2$term_name,gnlist_com_6v12v24h_dw2$p_value) 

# order data by p value
dotp_6h_up=dotp_6h_up[order(dotp_6h_up$gnlist_gnname_6h_up2.p_value),]
dotp_12h_up=dotp_12h_up[order(dotp_12h_up$gnlist_gnname_12h_up2.p_value),]
dotp_24h_up=dotp_24h_up[order(dotp_24h_up$gnlist_gnname_24h_up2.p_value),]
# order data by p value  down
dotp_6h_dw=dotp_6h_dw[order(dotp_6h_dw$gnlist_gnname_6h_dw2.p_value),]
dotp_12h_dw=dotp_12h_dw[order(dotp_12h_dw$gnlist_gnname_12h_dw2.p_value),]
dotp_24h_dw=dotp_24h_dw[order(dotp_24h_dw$gnlist_gnname_24h_dw2.p_value),]
# order common p value
dotp_6v12_up=dotp_6v12_up[order(dotp_6v12_up$gnlist_com_6v12h_up2.p_value),]
dotp_6v24_up=dotp_6v24_up[order(dotp_6v24_up$gnlist_com_6v24h_up2.p_value ),]
dotp_12v24_up=dotp_12v24_up[order(dotp_12v24_up$gnlist_com_12v24h_up2.p_value),]
dotp_6v12v24_up=dotp_6v12v24_up[order(dotp_6v12v24_up$gnlist_com_6v12v24h_up2.p_value),]
# order common dw
dotp_6v12_dw=dotp_6v12_dw[order(dotp_6v12_dw$gnlist_com_6v12h_dw2.p_value),]
dotp_6v24_dw=dotp_6v24_dw[order(dotp_6v24_dw$gnlist_com_6v24h_dw2.p_value),]
dotp_12v24_dw=dotp_12v24_dw[order(dotp_12v24_dw$gnlist_com_12v24h_dw2.p_value),]
dotp_6v12v24_dw=dotp_6v12v24_dw[order(dotp_6v12v24_dw$gnlist_com_6v12v24h_dw2.p_value),]


write.csv(dotp_6h_up, file = "genelist_pvalue_6_up_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12h_up, file = "genelist_pvalue_12_up_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_24h_up, file = "genelist_pvalue_24_up_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6h_dw, file = "genelist_pvalue_6_dw_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12h_dw, file = "genelist_pvalue_12_dw_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_24h_dw, file = "genelist_pvalue_24_dw_mi_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)

write.csv(dotp_6v12_up, file = "genelist_pvalue_6v12_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v24_up, file = "genelist_pvalue_6v24_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12v24_up, file = "genelist_pvalue_12v24_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v12v24_up, file = "genelist_pvalue_6v12v24_up_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v12_dw, file = "genelist_pvalue_6v12_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v24_dw, file = "genelist_pvalue_6v24_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_12v24_dw, file = "genelist_pvalue_12v24_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)
write.csv(dotp_6v12v24_dw, file = "genelist_pvalue_6v12v24_dw_ensg.csv", sep=",",eol = "\n", row.names=TRUE, col.names= TRUE)


# first transfer pvalue to numeric or do like up part
dotp_6h_dw$gnlist_gnname_6h_dw2.p_value= as.numeric(dotp_6h_dw$gnlist_gnname_6h_dw2.p_value)
dotp_12h_dw$gnlist_gnname_12h_dw2.p_value= as.numeric(dotp_12h_dw$gnlist_gnname_12h_dw2.p_value)
dotp_24h_dw$gnlist_gnname_24h_dw2.p_value= as.numeric(dotp_24h_dw$gnlist_gnname_24h_dw2.p_value)

# dot plot with order
ggplot(dotp_6h_up, aes(x = dotp_6h_up[,2], y = reorder(dotp_6h_up[,1],-as.numeric(dotp_6h_up[,2])))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#FFAAD4", stackratio = 1.5, dotsize = 0.5)+
  labs(title= "dotplot_6h_up_mi_ensg", x = "p value", y = NULL)
ggplot(dotp_12h_up, aes(x = dotp_12h_up[,2], y = reorder(dotp_12h_up[,1],-as.numeric(dotp_12h_up[,2]))))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#aafbff", stackratio = 1.5, dotsize = 0.5)+
  labs(title= "dotplot_12h_up_mi_ensg", x = "p value", y = NULL)
ggplot(dotp_24h_up, aes(x = dotp_24h_up[,2], y = reorder(dotp_24h_up[,1],-as.numeric(dotp_24h_up[,2]))))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#ffe1aa", stackratio = 1, dotsize = 0.5)+
  labs(title= "dotplot_24h_up_mi_ensg", x = "p value", y = NULL)

# dot plot with order down
ggplot(dotp_6h_dw, aes(x = dotp_6h_dw[,2], y = reorder(dotp_6h_dw[,1],-(dotp_6h_dw[,2])))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#FFAAD4", stackratio = 1, dotsize = 0.5)+
  labs(title= "dotplot_6h_dw_mi_ensg", x = "p value", y = NULL)
ggplot(dotp_12h_dw, aes(x = dotp_12h_dw[,2], y = reorder(dotp_12h_dw[,1],-dotp_12h_dw[,2]))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#aafbff", stackratio = 1, dotsize = 0.5)+
  labs(title= "dotplot_12h_dw_mi_ensg", x = "p value", y = NULL)
ggplot(dotp_24h_dw, aes(x = dotp_24h_dw[,2], y = reorder(dotp_24h_dw[,1],-dotp_24h_dw[,2]))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = "#ffe1aa", stackratio = 1, dotsize = 0.5)+
  labs(title= "dotplot_24h_dw_mi_ensg", x = "p value", y = NULL)

# save plot  
pdf("dotplot_6h_up_mi_ensg.pdf")
pdf("dotplot_12h_up_mi_ensg.pdf")
pdf("dotplot_24h_up_mi_ensg.pdf")
pdf("dotplot_6h_dw_mi_ensg.pdf")
pdf("dotplot_12h_dw_mi_ensg.pdf")
pdf("dotplot_24h_dw_mi_ensg.pdf")
dev.off()

#export genenames from every timepoint
miRNA_genenames= c()
miRNA_genenames= c(genenames_6h_up,genenames_12h_up,genenames_24h_up,genenames_6h_dw,genenames_12h_dw,genenames_24h_dw)
lapply(miRNA_genenames, function(x) write.table( data.frame(x), 'genenames_miRNA_ensg.csv'  , append= T, sep=',' ))

# export geneset term ID from every timepoint
gnlist_6h_up=as.data.frame(cbind(gnlist_gnname_6h_up2$term_ID, gnlist_gnname_6h_up2$term_name))
gnlist_12h_up=as.data.frame(cbind(gnlist_gnname_12h_up2$term_ID, gnlist_gnname_12h_up2$term_name))
gnlist_24h_up=as.data.frame(cbind(gnlist_gnname_24h_up2$term_ID, gnlist_gnname_24h_up2$term_name))
gnlist_6h_dw=as.data.frame(cbind(gnlist_gnname_6h_dw2$term_ID, gnlist_gnname_6h_dw2$term_name))
gnlist_12h_dw=as.data.frame(cbind(gnlist_gnname_12h_dw2$term_ID, gnlist_gnname_12h_dw2$term_name))
gnlist_24h_dw=as.data.frame(cbind(gnlist_gnname_24h_dw2$term_ID, gnlist_gnname_24h_dw2$term_name))
miRNA_gnlist= c(gnlist_6h_up,gnlist_12h_up,gnlist_24h_up,gnlist_6h_dw, gnlist_12h_dw, gnlist_24h_dw)
lapply(miRNA_gnlist, function(x) write.table( data.frame(x), 'geneset_miRNA_ensg.csv'  , append= T, sep=',' ))

#get mirbase gene name
mirbase_genenames_6h_up=biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name', "mirbase_id","hgnc_symbol"),
                                       filters = 'hgnc_symbol', 
                                       values = genenames_6h_up$hgnc_symbol,
                                       mart = mart)

mirbase_genenames_12h_up=biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name', "mirbase_id","hgnc_symbol"),
                                       filters = 'hgnc_symbol', 
                                       values = genenames_12h_up$hgnc_symbol,
                                       mart = mart)

mirbase_genenames_24h_up=biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name', "mirbase_id","hgnc_symbol"),
                                       filters = 'hgnc_symbol', 
                                       values = genenames_24h_up$hgnc_symbol,
                                       mart = mart)

mirbase_genenames_6h_dw=biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name', "mirbase_id","hgnc_symbol"),
                                       filters = 'hgnc_symbol', 
                                       values = genenames_6h_dw$hgnc_symbol,
                                       mart = mart)

mirbase_genenames_12h_dw=biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name', "mirbase_id","hgnc_symbol"),
                                        filters = 'hgnc_symbol', 
                                        values = genenames_12h_dw$hgnc_symbol,
                                        mart = mart)

mirbase_genenames_24h_dw=biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name', "mirbase_id","hgnc_symbol"),
                                        filters = 'hgnc_symbol', 
                                        values = genenames_24h_dw$hgnc_symbol,
                                        mart = mart)

#up
venn_list=list("mirbase"=mirbase_genenames_6h_up$hgnc_symbol,"gnname"=genenames_6h_up$hgnc_symbol)
venn.diagram(venn_list, filename = 'venn_6h_up.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, cat.col = c('red', 'blue'), 
             col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')
part_6h_up= get.venn.partitions(venn_list)
part_list_6h_up=list('set'=part$..set..,'count'= part$..count.., 'value'=part$..values..)

venn_list=list("mirbase"=mirbase_genenames_12h_up$hgnc_symbol,"gnname"=genenames_12h_up$hgnc_symbol)
venn.diagram(venn_list, filename = 'venn_12h_up.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, cat.col = c('red', 'blue'), 
             col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')
part_12h_up= get.venn.partitions(venn_list)
part_list_12h_up=list('set'=part$..set..,'count'= part$..count.., 'value'=part$..values..)


venn_list=list("mirbase"=mirbase_genenames_24h_up$hgnc_symbol,"gnname"=genenames_24h_up$hgnc_symbol)
venn.diagram(venn_list, filename = 'venn_24h_up.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, cat.col = c('red', 'blue'), 
             col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')
part_24h_up= get.venn.partitions(venn_list)
part_list_24h_up=list('set'=part$..set..,'count'= part$..count.., 'value'=part$..values..)

#dw
venn_list=list("mirbase"=mirbase_genenames_6h_dw$hgnc_symbol,"gnname"=genenames_6h_dw$hgnc_symbol)
venn.diagram(venn_list, filename = 'venn_6h_dw.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, cat.col = c('red', 'blue'), 
             col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')
part_6h_dw= get.venn.partitions(venn_list)
part_list_6h_dw=list('set'=part$..set..,'count'= part$..count.., 'value'=part$..values..)

venn_list=list("mirbase"=mirbase_genenames_12h_dw$hgnc_symbol,"gnname"=genenames_12h_dw$hgnc_symbol)
venn.diagram(venn_list, filename = 'venn_12h_dw.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, cat.col = c('red', 'blue'), 
             col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')
part_12h_dw= get.venn.partitions(venn_list)
part_list_12h_dw=list('set'=part$..set..,'count'= part$..count.., 'value'=part$..values..)

venn_list=list("mirbase"=mirbase_genenames_24h_dw$hgnc_symbol,"gnname"=genenames_24h_dw$hgnc_symbol)
venn.diagram(venn_list, filename = 'venn_24h_dw.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, cat.col = c('red', 'blue'), 
             col = c('red', 'blue'), cex = 1.5, fontfamily = 'serif', at.cex = 1.5, cat.fontfamily = 'serif')
part_24h_dw= get.venn.partitions(venn_list)
part_list_24h_dw=list('set'=part$..set..,'count'= part$..count.., 'value'=part$..values..)

# export geneset term ID from every timepoint
mirbase_gnlist_6h_up=as.data.frame(cbind(mirbase_genenames_6h_up$mirbase_id, mirbase_genenames_6h_up$hgnc_symbol))
mirbase_gnlist_12h_up=as.data.frame(cbind(mirbase_genenames_12h_up$mirbase_id, mirbase_genenames_12h_up$hgnc_symbol))
mirbase_gnlist_24h_up=as.data.frame(cbind(mirbase_genenames_24h_up$mirbase_id, mirbase_genenames_24h_up$hgnc_symbol))
mirbase_gnlist_6h_dw=as.data.frame(cbind(mirbase_genenames_6h_dw$mirbase_id, mirbase_genenames_6h_dw$hgnc_symbol))
mirbase_gnlist_12h_dw=as.data.frame(cbind(mirbase_genenames_12h_dw$mirbase_id, mirbase_genenames_12h_dw$hgnc_symbol))
mirbase_gnlist_24h_dw=as.data.frame(cbind(mirbase_genenames_24h_dw$mirbase_id, mirbase_genenames_24h_dw$hgnc_symbol))
miRNA_mirbase_gnlist= c(mirbase_gnlist_6h_up,mirbase_gnlist_12h_up,mirbase_gnlist_24h_up,mirbase_gnlist_6h_dw, mirbase_gnlist_12h_dw, mirbase_gnlist_24h_dw)
lapply(miRNA_mirbase_gnlist, function(x) write.table( data.frame(x), 'mirbase_gnlist_miRNA.csv'  , append= T, sep=',' ))


#compare human gene symbol to mirna
setwd("/home/yuanpeng/yuan/R sarsCoV2/mirna_sequen/predict")
targetscan_dock <- read.table(file = "TargetScan-dock1-predict-mirna.csv", header = TRUE, sep = ",")
targetscan_RAB11FIP1 <- read.table(file = "TargetScan-RAB11FIP1-predict-mirna.csv", header = TRUE, sep = "\t")

setwd("/home/yuanpeng/yuan/R sarsCoV2/mirna_sequen")
mirabase_gnlist <- read.table(file = "mirbase_gnlist_miRNA.csv", header = TRUE, sep = ",")

#>>
mirabase_list=list("mirbase_6h_up"=mirabase_gnlist$mirbase_6h_up,"hgnc_6h_up"=mirabase_gnlist$hgnc_6h_up,"mirbase_12h_up"=mirabase_gnlist$mirbase_12h_up,"hgnc_12h_up"=mirabase_gnlist$hgnc_12h_up,
                   "mirbase_24h_up"=mirabase_gnlist$mirbase_24h_up,"hgnc_24h_up"=mirabase_gnlist$hgnc_24h_up,"mirbase_6h_dw"=mirabase_gnlist$mirbase_6h_dw,"hgnc_6h_dw"=mirabase_gnlist$hgnc_6h_dw,
                   "mirbase_12h_dw"=mirabase_gnlist$mirbase_12h_dw,"hgnc_12h_dw"=mirabase_gnlist$hgnc_12h_dw,"mirbase_24h_dw"=mirabase_gnlist$mirbase_24h_dw,"hgnc_24h_dw"=mirabase_gnlist$hgnc_24h_dw)
sapply(mirabase_list, function(y) "" %in% y)
mirabase_list$mirbase_6h_up = mirabase_list$mirbase_6h_up [mirabase_list$mirbase_6h_up  != ""]
mirabase_list$hgnc_6h_up = mirabase_list$hgnc_6h_up [mirabase_list$hgnc_6h_up  != ""]
mirabase_list$mirbase_12h_up = mirabase_list$mirbase_12h_up [mirabase_list$mirbase_12h_up  != ""]
mirabase_list$hgnc_12h_up = mirabase_list$hgnc_12h_up [mirabase_list$hgnc_12h_up  != ""]
mirabase_list$mirbase_24h_up = mirabase_list$mirbase_24h_up [mirabase_list$mirbase_24h_up  != ""]
mirabase_list$hgnc_24h_up = mirabase_list$hgnc_24h_up [mirabase_list$hgnc_24h_up  != ""]

mirabase_list$mirbase_6h_dw = mirabase_list$mirbase_6h_dw [mirabase_list$mirbase_6h_dw  != ""]
mirabase_list$hgnc_6h_dw = mirabase_list$hgnc_6h_dw [mirabase_list$hgnc_6h_dw  != ""]
mirabase_list$mirbase_12h_dw = mirabase_list$mirbase_12h_dw [mirabase_list$mirbase_12h_  != ""]
mirabase_list$hgnc_12h_dw = mirabase_list$hgnc_12h_dw [mirabase_list$hgnc_12h_dw  != ""]
mirabase_list$mirbase_24h_dw = mirabase_list$mirbase_24h_dw [mirabase_list$mirbase_24h_dw  != ""]
mirabase_list$hgnc_24h_dw = mirabase_list$hgnc_24h_dw [mirabase_list$hgnc_24h_dw  != ""]
#>>

intersect(targetscan_dock$miRNA,mirabase_gnlist$mirbase_6h_up)
intersect(targetscan_dock$miRNA,mirabase_gnlist$mirbase_12h_up)
intersect(targetscan_dock$miRNA,mirabase_gnlist$mirbase_24h_up)
intersect(targetscan_dock$miRNA,mirabase_gnlist$mirbase_6h_dw)
intersect(targetscan_dock$miRNA,mirabase_gnlist$mirbase_12h_dw)
intersect(targetscan_dock$miRNA,mirabase_gnlist$mirbase_24h_dw)

intersect(targetscan_RAB11FIP1$miRNA,mirabase_gnlist$mirbase_6h_up)
intersect(targetscan_RAB11FIP1$miRNA,mirabase_gnlist$mirbase_12h_up)
intersect(targetscan_RAB11FIP1$miRNA,mirabase_gnlist$mirbase_24h_up)
intersect(targetscan_RAB11FIP1$miRNA,mirabase_gnlist$mirbase_6h_dw)
intersect(targetscan_RAB11FIP1$miRNA,mirabase_gnlist$mirbase_12h_dw)
intersect(targetscan_RAB11FIP1$miRNA,mirabase_gnlist$mirbase_24h_dw)

