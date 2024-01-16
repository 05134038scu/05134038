# gn redundandecy
# genesetRedundancy
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(gurobi)
library(slam)
library(ontologyLIP)

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
dwgn_mrna_gnset_100top <- read.csv(file = "dwgn_mrna_gnset_100top.csv", header = TRUE, sep = ",")
dwgn_cir_gnset_100top<- read.csv(file = "dwgn_cir_gnset_100top.csv", header = TRUE, sep = ",")
upgn_cir_gnset_100top<- read.csv(file = "upgn_cir_gnset_100top.csv", header = TRUE, sep = ",")


setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/ontol")

tem_data<- data.frame("term_ID" = as.character(dwgn_mrna_gnset_100top$term_ID),
                      "p_value" = dwgn_mrna_gnset_100top$p_value,
                      "intersection" = as.character(dwgn_mrna_gnset_100top$intersection))

tem_data$term_ID <- as.character(tem_data$term_ID)
tem_data$intersection <- as.character(tem_data$intersection)
# first graph 
geneFilteringProfile(tem_data,plotName = "geneFiltering_dw_mrna_gnset_100top")
# second graph 
genesetRedundancyProfile(tem_data, plotName = "genesetRedundancy_dw_mrna_gnset_100top")
results <- removeRedundant(tem_data, jacCutoff = 0.3)

Pos <- sapply(results$term_ID, function(x) which(x == dwgn_mrna_gnset_100top$term_ID))
results_jacCut <- dwgn_mrna_gnset_100top[unlist(Pos),]
write.table(results_jacCut, "ontologyLIPresult_dw_mrna_gnset_100top-0.3.csv", sep = ",")

# dw cir   0.5 89% 其他都100% 所以沒有redun
tem_data<- data.frame("term_ID" = as.character(dwgn_cir_gnset_100top$term_ID),
                      "p_value" = dwgn_cir_gnset_100top$p_value,
                      "intersection" = as.character(dwgn_cir_gnset_100top$intersection))

tem_data$term_ID <- as.character(tem_data$term_ID)
tem_data$intersection <- as.character(tem_data$intersection)
# first graph 
geneFilteringProfile(tem_data,plotName = "geneFiltering_dw_cir_gnset_100top")
# second graph 
genesetRedundancyProfile(tem_data, plotName = "genesetRedundancy_dw_cir_gnset_100top")
results <- removeRedundant(tem_data, jacCutoff = 0.3)

Pos <- sapply(results$term_ID, function(x) which(x == dwgn_cir_gnset_100top$term_ID))
results_jacCut <- dwgn_cir_gnset_100top[unlist(Pos),]
write.table(results_jacCut, "ontologyLIPresult_dw_cir_gnset_100top-0.3.csv", sep = ",")

# cir up  都100% 所以沒有redun
tem_data<- data.frame("term_ID" = as.character(upgn_cir_gnset_100top$term_ID),
                      "p_value" = upgn_cir_gnset_100top$p_value,
                      "intersection" = as.character(upgn_cir_gnset_100top$intersection))

tem_data$term_ID <- as.character(tem_data$term_ID)
tem_data$intersection <- as.character(tem_data$intersection)
# first graph 
geneFilteringProfile(tem_data,plotName = "geneFiltering_up_cir_gnset_100top")
# second graph 
genesetRedundancyProfile(tem_data, plotName = "genesetRedundancy_up_cir_gnset_100top")
results <- removeRedundant(tem_data, jacCutoff = 0.3)

Pos <- sapply(results$term_ID, function(x) which(x == upgn_cir_gnset_100top$term_ID))
results_jacCut <- upgn_cir_gnset_100top[unlist(Pos),]
write.table(results_jacCut, "ontologyLIPresult_up_cir_gnset_100top-0.3.csv", sep = ",")

#作 reannotate>>>>>>>>>>>>>>>>>>>>>>
# grep GO from gmt file from gprofiler website
#cant save so add , stringsAsFactors=FALSE
#save txt and change to .txt > .gmt 
#col name and row name should be false


setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reannotate")
gprofiler_full_hsapiens.name= read.csv(file = "gprofiler_full_hsapiens.name.gmt", header = F, sep = "\t")

#grep
gprofiler_hsapiens_GO=gprofiler_full_hsapiens.name[grep("GO",gprofiler_full_hsapiens.name$V1),]
newGO=c()

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
upgn_mrna_gnset_100top2 = read.table(file = "upgn_mrna_gnset_100top.csv", header = TRUE, sep = ",")
dwgn_mrna_gnset_100top2 = read.table(file = "dwgn_mrna_gnset_100top.csv", header = TRUE, sep = ",")
upgn_cir_gnset_100top2 = read.table(file = "upgn_cir_gnset_100top.csv", header = TRUE, sep = ",")
dwgn_cir_gnset_100top2 = read.table(file = "dwgn_cir_gnset_100top.csv", header = TRUE, sep = ",")

totalgnset=c()
totalgnset= bind_rows(upgn_mrna_gnset_100top2,dwgn_mrna_gnset_100top2,
                    upgn_cir_gnset_100top2,dwgn_cir_gnset_100top2)
totalgnset=unique(totalgnset$term_ID)

newGO=c()
newGO= gprofiler_hsapiens_GO[gprofiler_hsapiens_GO$V1 %in% totalgnset,]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top100_reanno")
write.table(newGO,"gprofiler_hsapiens_newGO.csv", sep = "\t",col.names=FALSE, row.names = FALSE)


# test  token "gp__bdPA_Arny_ntk"

test=gprofiler2::gost(c("MAD1L1","POLDIP2"), correction_method="fdr",
                                        sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__bdPA_Arny_ntk" )


#reannoate top 100 gn gsea
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
upgn_mrna_ID= read.table(file = "upgn_mrna_gnname.csv", header = TRUE, sep = ",")
dwgn_mrna_ID= read.table(file = "dwgn_mrna_gnname.csv", header = TRUE, sep = ",")
upgn_cir_ID= read.table(file = "upgn_cir_gnname.csv", header = TRUE, sep = ",")
dwgn_cir_ID= read.table(file = "dwgn_cir_gnname.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
freq_mrna_100top= read.table(file = "freq_mrna_100top.csv", header = TRUE, sep = ",")
freq_cir_100top=read.table(file = "freq_cir_100top.csv", header = TRUE, sep = ",")

#拿已經有的 up dw gn去根 top100 比 得到其中的 up dw gn在作 GSEA
freq_mrna_up_100top= freq_mrna_100top[freq_mrna_100top$Var1 %in% upgn_mrna_ID$hgnc_symbol, ]
freq_mrna_dw_100top= freq_mrna_100top[freq_mrna_100top$Var1 %in% dwgn_mrna_ID$hgnc_symbol, ]
freq_cir_up_100top= freq_cir_100top[freq_cir_100top$Var1 %in% upgn_cir_ID$hgnc_symbol, ]
freq_cir_dw_100top= freq_cir_100top[freq_cir_100top$Var1 %in% dwgn_cir_ID$hgnc_symbol, ]


upgn_mrna_gnset_100top3 =gprofiler2::gost(freq_mrna_up_100top$Var1, correction_method="fdr",
                                            sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__bdPA_Arny_ntk" )
dwgn_mrna_gnset_100top3 =gprofiler2::gost(freq_mrna_dw_100top$Var1, correction_method="fdr",
                                          sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__bdPA_Arny_ntk" )
upgn_cir_gnset_100top3 =gprofiler2::gost(freq_cir_up_100top$Var1, correction_method="fdr",
                                          sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__bdPA_Arny_ntk" )
dwgn_cir_gnset_100top3 =gprofiler2::gost(freq_cir_dw_100top$Var1, correction_method="fdr",
                                          sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__bdPA_Arny_ntk" )

upgn_mrna_gnset_100top3=upgn_mrna_gnset_100top3$result
dwgn_mrna_gnset_100top3=dwgn_mrna_gnset_100top3$result
upgn_cir_gnset_100top3=upgn_cir_gnset_100top3$result
dwgn_cir_gnset_100top3=dwgn_cir_gnset_100top3$result



# 1116 2023 在作一次
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reannotate")
gprofiler_full_hsapiens.name= read.csv(file = "gprofiler_full_hsapiens.name.gmt", header = F, sep = "\t")

#grep
gprofiler_hsapiens_GO=gprofiler_full_hsapiens.name[grep("GO",gprofiler_full_hsapiens.name$V1),]


setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
upgn_mrna_gnset_100top2 = read.table(file = "upgn_mrna_gnset_100top.csv", header = TRUE, sep = ",")
upgn_cir_gnset_100top2 = read.table(file = "upgn_cir_gnset_100top.csv", header = TRUE, sep = ",")
dwgn_cir_gnset_100top2 = read.table(file = "dwgn_cir_gnset_100top.csv", header = TRUE, sep = ",")
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/ontol")
dwgn_mrna_gnset_100top2 = read.table(file = "ontologyLIPresult_dw_mrna_gnset_100top-0.3.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top100_reanno")

unique_total_gnlist_updw= c(upgn_mrna_gnset_100top2$term_ID,dwgn_mrna_gnset_100top2 $term_ID,
                            upgn_cir_gnset_100top2$term_ID,dwgn_cir_gnset_100top2$term_ID)
unique_total_gnlist_updw= unique(unique_total_gnlist_updw)

newGO=c()
newGO= gprofiler_hsapiens_GO[gprofiler_hsapiens_GO$V1 %in% unique_total_gnlist_updw,]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top100_reanno")
write.table(newGO,"gprofiler_hsapiens_newGO222.csv", sep = "\t",col.names=FALSE, row.names = FALSE)


#  gp__JNp9_kc8g_47s
#reannoate top 100 gn gsea
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
upgn_mrna_ID= read.table(file = "upgn_mrna_gnname.csv", header = TRUE, sep = ",")
dwgn_mrna_ID= read.table(file = "dwgn_mrna_gnname.csv", header = TRUE, sep = ",")
upgn_cir_ID= read.table(file = "upgn_cir_gnname.csv", header = TRUE, sep = ",")
dwgn_cir_ID= read.table(file = "dwgn_cir_gnname.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
freq_mrna_100top= read.table(file = "freq_mrna_100top.csv", header = TRUE, sep = ",")
freq_cir_100top=read.table(file = "freq_cir_100top.csv", header = TRUE, sep = ",")

#拿已經有的 up dw gn去根 top100 比 得到其中的 up dw gn在作 GSEA
freq_mrna_up_100top= freq_mrna_100top[freq_mrna_100top$Var1 %in% upgn_mrna_ID$hgnc_symbol, ]
freq_mrna_dw_100top= freq_mrna_100top[freq_mrna_100top$Var1 %in% dwgn_mrna_ID$hgnc_symbol, ]
freq_cir_up_100top= freq_cir_100top[freq_cir_100top$Var1 %in% upgn_cir_ID$hgnc_symbol, ]
freq_cir_dw_100top= freq_cir_100top[freq_cir_100top$Var1 %in% dwgn_cir_ID$hgnc_symbol, ]


upgn_mrna_gnset_100top3 =gprofiler2::gost(freq_mrna_up_100top$Var1, correction_method="fdr",
                                          sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__JNp9_kc8g_47s" )
dwgn_mrna_gnset_100top3 =gprofiler2::gost(freq_mrna_dw_100top$Var1, correction_method="fdr",
                                          sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__JNp9_kc8g_47s" )
upgn_cir_gnset_100top3 =gprofiler2::gost(freq_cir_up_100top$Var1, correction_method="fdr",
                                         sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__JNp9_kc8g_47s" )
dwgn_cir_gnset_100top3 =gprofiler2::gost(freq_cir_dw_100top$Var1, correction_method="fdr",
                                         sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism ="gp__JNp9_kc8g_47s" )

upgn_mrna_gnset_100top3=upgn_mrna_gnset_100top3$result
dwgn_mrna_gnset_100top3=dwgn_mrna_gnset_100top3$result
upgn_cir_gnset_100top3=upgn_cir_gnset_100top3$result
dwgn_cir_gnset_100top3=dwgn_cir_gnset_100top3$result  #找不到gnset

upgn_mrna_gnset_100top4= upgn_mrna_gnset_100top3[,c(3:6,9,11,16)]
dwgn_mrna_gnset_100top4= dwgn_mrna_gnset_100top3[,c(3:6,9,11,16)]
upgn_cir_gnset_100top4= upgn_cir_gnset_100top3[,c(3:6,9,11,16)]

dwgn_mrna_gnset_100top4=dwgn_mrna_gnset_100top4[order(as.numeric(dwgn_mrna_gnset_100top4$p_value)),]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/ontol")
write.csv(upgn_mrna_gnset_100top4, file = "upgn_mrna_gnset_100top_reanno.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_mrna_gnset_100top4, file = "dwgn_mrna_gnset_100top_reanno.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(upgn_cir_gnset_100top4, file = "upgn_cir_gnset_100top_reanno.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# reanno gorup
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")

annodata= read.csv(file = "AllAnnoData3.csv", header = TRUE, sep = "\t",quote = "", stringsAsFactors=FALSE )

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/reanno_group")

upgn_mrna_gnset_100top_reanno= read.table(file = "upgn_mrna_gnset_100top_reanno.csv", header = TRUE, sep = ",")

upgn_mrna_gnset_100top_group= merge(x=upgn_mrna_gnset_100top4, y=annodata, by.x=c("term_name"), by.y = c("term_name"))
dwgn_mrna_gnset_100top_group= merge(x=dwgn_mrna_gnset_100top4, y=annodata, by.x=c("term_name"), by.y = c("term_name"))
upgn_cir_gnset_100top_group= merge(x=upgn_cir_gnset_100top4, y=annodata, by.x=c("term_name"), by.y = c("term_name"))

# not in

upgn_mrna_gnset_100top_group_not= upgn_mrna_gnset_100top4[!upgn_mrna_gnset_100top4$term_name %in% annodata$term_name,]
dwgn_mrna_gnset_100top_group_not= dwgn_mrna_gnset_100top4[!dwgn_mrna_gnset_100top4$term_name %in% annodata$term_name,]
upgn_cir_gnset_100top_group_not= upgn_cir_gnset_100top4[!upgn_mrna_gnset_100top4$term_name %in% annodata$term_name,]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/reanno_group/notgroup")
write.csv(upgn_mrna_gnset_100top_group_not, file = "upgn_mrna_gnset_100top_group_not.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_mrna_gnset_100top_group_not, file = "dwgn_mrna_gnset_100top_group_not.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(upgn_cir_gnset_100top_group_not, file = "upgn_cir_gnset_100top_group_not.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# 把不再annodata的 term 加進table後 在merge一次
upgn_mrna_gnset_100top_group=c()
dwgn_mrna_gnset_100top_group=c()
upgn_cir_gnset_100top_group=c()

upgn_mrna_gnset_100top_group= merge(x=upgn_mrna_gnset_100top4, y=annodata, by.x=c("term_name"), by.y = c("term_name"))
dwgn_mrna_gnset_100top_group= merge(x=dwgn_mrna_gnset_100top4, y=annodata, by.x=c("term_name"), by.y = c("term_name"))
upgn_cir_gnset_100top_group= merge(x=upgn_cir_gnset_100top4, y=annodata, by.x=c("term_name"), by.y = c("term_name"))

dwgn_mrna_gnset_100top_group=dwgn_mrna_gnset_100top_group[!duplicated(dwgn_mrna_gnset_100top_group$term_name),]
upgn_cir_gnset_100top_group=upgn_cir_gnset_100top_group[!duplicated(upgn_cir_gnset_100top_group$term_name),]

dwgn_mrna_gnset_100top_group=dwgn_mrna_gnset_100top_group[order(as.numeric(dwgn_mrna_gnset_100top_group$p_value)),]
upgn_cir_gnset_100top_group=upgn_cir_gnset_100top_group[order(as.numeric(upgn_cir_gnset_100top_group$p_value)),]


setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/reanno_group")
write.csv(upgn_mrna_gnset_100top_group, file = "upgn_mrna_gnset_100top_reanno_group.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_mrna_gnset_100top_group, file = "dwgn_mrna_gnset_100top_reanno_group.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(upgn_cir_gnset_100top_group, file = "upgn_cir_gnset_100top_reanno_group.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


#看原本900組的 rna 根 top 100 node 的共同gnset有幾個
# make geneset up-up list and dw-dw list
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/ontol_up")
upgn_mrna_gnset2= read.table(file = "upgn_mrna_gnset_reanno.csv", header = TRUE, sep = ",")
upgn_cir_gnset2= read.table(file = "upgn_cir_gnset_reanno.csv", header = TRUE, sep = ",")
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/ontol_dw")
dwgn_mrna_gnset2= read.table(file = "dwgn_mrna_gnset_reanno.csv", header = TRUE, sep = ",")
dwgn_cir_gnset2= read.table(file = "dwgn_cir_gnset(ori).csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/reanno_group")
upgn_mrna_gnset_100top_group= read.table(file = "upgn_mrna_gnset_100top_reanno_group.csv", header = TRUE, sep = ",")
dwgn_mrna_gnset_100top_group= read.table(file = "dwgn_mrna_gnset_100top_reanno_group.csv", header = TRUE, sep = ",")
upgn_cir_gnset_100top_group= read.table(file = "upgn_cir_gnset_100top_reanno_group.csv", header = TRUE, sep = ",")


# 取得 common gnset 1116 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/common/common2")
up_up_mrna= list("up_mRNA_reanno"=upgn_mrna_gnset2$term_ID,"up_mRNA_top100"=upgn_mrna_gnset_100top_group$term_id)
up_up_cir= list("up_cir_reanno"=upgn_cir_gnset2$term_ID,"up_cir_top100"=upgn_cir_gnset_100top_group$term_id)

dw_dw_mrna= list("dw_mRNA_reanno"=dwgn_mrna_gnset2$term_ID,"dw_mRNA_top100"=dwgn_mrna_gnset_100top_group$term_id)

# venndiagram
venn.diagram(up_up_mrna, filename = 'up_up_mrna_gnset.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'),col = c('red', 'blue'), cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif',margin=0.05,cat.pos=0)
venn.diagram(up_up_cir, filename = 'up_up_cir_gnset.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'),col = c('red', 'blue'), cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif',margin=0.05,cat.pos=0)
venn.diagram(dw_dw_mrna, filename = 'dw_dw_mrna_gnset.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'),col = c('red', 'blue'), cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif',margin=0.05,cat.pos=0) #,cat.dist=0.03

# 取得 common gnset 1102 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

up_up_cir_gnset= merge(x=upgn_cir_gnset2, y=upgn_cir_gnset_100top_group, by.x=c("term_ID"), by.y = c("term_id"))
dw_dw_mrna_gnset=  merge(x=dwgn_mrna_gnset2, y=dwgn_mrna_gnset_100top_group, by.x=c("term_ID"), by.y = c("term_id"))



#在首動 作表格  
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/common/common2")

write.csv(up_up_cir_gnset, file = "up_up_cir_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dw_dw_mrna_gnset, file = "dw_dw_mrna_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)



# 加上顏色 作長條圖
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")
color_code= read.csv(file = "color_code.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/common/common2")
dw_dw_mrna_gnset= read.csv(file = "dw_dw_mrna_gnset.csv", header = TRUE, sep = ",")

# dwgn_mrna_gnset3 又部能讀檔案 去用python作   mrna序列且mi-mrna結果.ipynb
# 發現用read csv就可以讀了
dw_dw_mrna_gnset_group= dw_dw_mrna_gnset %>%
  group_by(group) %>%
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

dw_dw_mrna_gnset_group2= merge(x=dw_dw_mrna_gnset_group, y=color_code, by.x=c("group"), by.y = c("group"))

category_counts <- table(dw_dw_mrna_gnset$group)
category_counts <- sort(category_counts, decreasing = TRUE)
category_counts <- data.frame(category_counts)
s <- sum(category_counts$Freq)
category_counts$Freq <- category_counts$Freq / sum(category_counts$Freq)*100

category_counts2= merge(x=category_counts, y=color_code, by.x=c("Var1"), by.y = c("group"))
category_counts2=category_counts2[order(category_counts2$Freq, decreasing = TRUE),]

ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=round(as.numeric(Freq), digits = 1)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip(ylim=c(0,65))  + #把圖轉至
  labs( x = "", y = "",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic(base_size = 15)  + #把所有原本字體同時設定
  theme()+
  theme(legend.position = "none") #要有這個把legend刪掉

# 取得 common gnset 1116 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


