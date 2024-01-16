# 先作 onotology lip
# genesetRedundancy
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(gurobi)
library(slam)
library(ontologyLIP)

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
upgn_mrna_gnset= read.csv(file = "upgn_mrna_gnset.csv", header = T, sep = ",")
dwgn_mrna_gnset= read.csv(file = "dwgn_mrna_gnset.csv", header = T, sep = ",")
upgn_cir_gnset= read.csv(file = "upgn_cir_gnset.csv", header = T, sep = ",")
dwgn_cir_gnset= read.csv(file = "dwgn_cir_gnset.csv", header = T, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")

# mrna
tem_data<- data.frame("term_ID" = as.character(upgn_mrna_gnset$term_ID),
                      "p_value" = upgn_mrna_gnset$p_value,
                      "intersection" = as.character(upgn_mrna_gnset$intersection))

tem_data$term_ID <- as.character(tem_data$term_ID)
tem_data$intersection <- as.character(tem_data$intersection)

# first graph for cutoff
geneFilteringProfile(tem_data,plotName = "geneFiltering_up_mrna")

# second graph to see how the terms get reduced
genesetRedundancyProfile(tem_data, plotName = "genesetRedundancy_up_mrna")

# write the cutoff here keep 0.95 remain at least
results <- removeRedundant(tem_data, jacCutoff = 0.2)

Pos <- sapply(results$term_ID, function(x) which(x == upgn_mrna_gnset$term_ID))
results_jacCut <- upgn_mrna_gnset[unlist(Pos),]

write.table(results_jacCut, "ontologyLIPresult_up_mrna-0.2.csv", sep = ",")


tem_data<- data.frame("term_ID" = as.character(dwgn_mrna_gnset$term_ID),
                      "p_value" = dwgn_mrna_gnset$p_value,
                      "intersection" = as.character(dwgn_mrna_gnset$intersection))

tem_data$term_ID <- as.character(tem_data$term_ID)
tem_data$intersection <- as.character(tem_data$intersection)
# first graph 
geneFilteringProfile(tem_data,plotName = "geneFiltering_dw_mrna")
# second graph 
genesetRedundancyProfile(tem_data, plotName = "genesetRedundancy_dw_mrna")
results <- removeRedundant(tem_data, jacCutoff = 0.2)

Pos <- sapply(results$term_ID, function(x) which(x == dwgn_mrna_gnset$term_ID))
results_jacCut <- dwgn_mrna_gnset[unlist(Pos),]
write.table(results_jacCut, "ontologyLIPresult_dw_mrna-0.2.csv", sep = ",")

# cir
tem_data<- data.frame("term_ID" = as.character(upgn_cir_gnset$term_ID),
                      "p_value" = upgn_cir_gnset$p_value,
                      "intersection" = as.character(upgn_cir_gnset$intersection))

tem_data$term_ID <- as.character(tem_data$term_ID)
tem_data$intersection <- as.character(tem_data$intersection)
# first graph 
geneFilteringProfile(tem_data,plotName = "geneFiltering_up_cir")
# second graph 
genesetRedundancyProfile(tem_data, plotName = "genesetRedundancy_up_cir")
results <- removeRedundant(tem_data, jacCutoff = 0.2)

Pos <- sapply(results$term_ID, function(x) which(x == upgn_cir_gnset$term_ID))
results_jacCut <- upgn_cir_gnset[unlist(Pos),]
write.table(results_jacCut, "ontologyLIPresult_up_cir-0.2.csv", sep = ",")


tem_data<- data.frame("term_ID" = as.character(dwgn_cir_gnset$term_ID),
                      "p_value" = dwgn_cir_gnset$p_value,
                      "intersection" = as.character(dwgn_cir_gnset$intersection))

tem_data$term_ID <- as.character(tem_data$term_ID)
tem_data$intersection <- as.character(tem_data$intersection)
# first graph 
geneFilteringProfile(tem_data,plotName = "geneFiltering_dw_cir")
# second graph 
genesetRedundancyProfile(tem_data, plotName = "genesetRedundancy_dw_cir")
results <- removeRedundant(tem_data, jacCutoff = 0.2)

Pos <- sapply(results$term_ID, function(x) which(x == dwgn_cir_gnset$term_ID))
results_jacCut <- dwgn_cir_gnset[unlist(Pos),]
write.table(results_jacCut, "ontologyLIPresult_dw_cir-0.2.csv", sep = ",")



#先 取得所有GO id from ontol result 作出新的 gmt file

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reannotate")
gprofiler_full_hsapiens.name= read.csv(file = "gprofiler_full_hsapiens.name.gmt", header = F, sep = "\t")
gprofiler_hsapiens_GO=gprofiler_full_hsapiens.name[grep("GO",gprofiler_full_hsapiens.name$V1),]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/ontol_up")
upgn_mrna_gnset= read.csv(file = "ontologyLIPresult_up_mrna-0.2.csv", header = T, sep = ",")
upgn_cir_gnset= read.csv(file = "ontologyLIPresult_up_cir-0.2.csv", header = T, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/ontol_dw")
dwgn_mrna_gnset= read.csv(file = "ontologyLIPresult_dw_mrna-0.2.csv", header = T, sep = ",")
dwgn_cir_gnset= read.csv(file = "ontologyLIPresult_dw_cir-0.2.csv", header = T, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reannotate")

unique_total_gnlist_updw= c(upgn_mrna_gnset$term_ID,dwgn_mrna_gnset$term_ID,upgn_cir_gnset$term_ID,dwgn_cir_gnset$term_ID)
unique_total_gnlist_updw= unique(unique_total_gnlist_updw)

newGO= gprofiler_hsapiens_GO[gprofiler_hsapiens_GO$V1 %in% unique_total_gnlist_updw,]

gprofiler_hsapiens_GO=data.frame(lapply(gprofiler_hsapiens_GO, as.character), stringsAsFactors=FALSE)

write.table(newGO,"gprofiler_hsapiens_newGO2.txt", sep = "\t",col.names=FALSE, row.names = FALSE)

#到linux 刪除 " 不然根本不能用 雖然有token 但是無法作GSEA
# sed -i 's/"//g' gprofiler_hsapiens_newGO3.txt


#手動改檔名 gmt 把這檔案gmt 丟回gprofiler 網頁上傳 custom gmt
# 丟到gmt helper vaildate 得到 token 在來跑 gprofiler
# g:Profiler: gp__CPLN_wZuy_rD4

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
upgn_mrna_gnname= read.table(file = "upgn_mrna_gnname.csv", header = TRUE, sep = ",")
dwgn_mrna_gnname=read.table(file = "dwgn_mrna_gnname.csv", header = TRUE, sep = ",")
upgn_cir_gnname= read.table(file = "upgn_cir_gnname.csv", header = TRUE, sep = ",")
dwgn_cir_gnname= read.table(file = "dwgn_cir_gnname.csv", header = TRUE, sep = ",")

df=gprofiler2::gost(query = c("KDM1A","AURKC","MRE11"), correction_method="fdr", 
                    sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism = "gp__CPLN_wZuy_rD4")

upgn_mrna_gnset=gprofiler2::gost(upgn_mrna_gnname$hgnc_symbol, correction_method="fdr", 
                                 sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism = "gp__CPLN_wZuy_rD4")
dwgn_mrna_gnset=gprofiler2::gost(dwgn_mrna_gnname$hgnc_symbol, correction_method="fdr", 
                                 sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism = "gp__CPLN_wZuy_rD4")
upgn_cir_gnset=gprofiler2::gost(upgn_cir_gnname$hgnc_symbol, correction_method="fdr", 
                                 sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism = "gp__CPLN_wZuy_rD4")
dwgn_cir_gnset=gprofiler2::gost(dwgn_cir_gnname$hgnc_symbol, correction_method="fdr", 
                                 sources =c("GO:BP","GO:MF"),evcodes = TRUE,organism = "gp__CPLN_wZuy_rD4")

upgn_mrna_gnset=upgn_mrna_gnset$result
dwgn_mrna_gnset=dwgn_mrna_gnset$result
upgn_cir_gnset=upgn_cir_gnset$result
dwgn_cir_gnset=dwgn_cir_gnset$result

# dwgn_cir_gnset=dwgn_cir_gnset$result  only have 10 gnset so here is no result

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")

upgn_mrna_gnset2<- data.frame()
for (i in 1:length(upgn_mrna_gnset$term_size)){
  if(upgn_mrna_gnset$term_size[i]<1000 && upgn_mrna_gnset$term_size[i]>10){
    print(i)
    print(upgn_mrna_gnset$term_size[i])
    upgn_mrna_gnset2=rbind(upgn_mrna_gnset2,c(upgn_mrna_gnset$term_size[i],
                                              upgn_mrna_gnset$p_value[i],upgn_mrna_gnset$term_id[i]
                                                    ,upgn_mrna_gnset$term_name[i],upgn_mrna_gnset$query_size[i]
                                                    ,upgn_mrna_gnset$intersection_size[i],upgn_mrna_gnset$intersection[i]))
  }
}

colnames(upgn_mrna_gnset2 )=c('term_size',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
upgn_mrna_gnset2=upgn_mrna_gnset2[order(as.numeric(upgn_mrna_gnset2$p_value)),]

write.csv(upgn_mrna_gnset2, file = "upgn_mrna_gnset_reanno.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

dwgn_mrna_gnset2<- data.frame()
for (i in 1:length(dwgn_mrna_gnset$term_size)){
  if(dwgn_mrna_gnset$term_size[i]<1000 && dwgn_mrna_gnset$term_size[i]>10){
    print(i)
    print(dwgn_mrna_gnset$term_size[i])
    dwgn_mrna_gnset2=rbind(dwgn_mrna_gnset2,c(dwgn_mrna_gnset$term_size[i],
                                              dwgn_mrna_gnset$p_value[i],dwgn_mrna_gnset$term_id[i]
                                              ,dwgn_mrna_gnset$term_name[i],dwgn_mrna_gnset$query_size[i]
                                              ,dwgn_mrna_gnset$intersection_size[i],dwgn_mrna_gnset$intersection[i]))
  }
}

colnames(dwgn_mrna_gnset2 )=c('term_size',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
dwgn_mrna_gnset2=dwgn_mrna_gnset2[order(as.numeric(dwgn_mrna_gnset2$p_value)),]

write.csv(dwgn_mrna_gnset2, file = "dwgn_mrna_gnset_reanno.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

upgn_cir_gnset2<- data.frame()
for (i in 1:length(upgn_cir_gnset$term_size)){
  if(upgn_cir_gnset$term_size[i]<1000 && upgn_cir_gnset$term_size[i]>10){
    print(i)
    print(upgn_cir_gnset$term_size[i])
    upgn_cir_gnset2=rbind(upgn_cir_gnset2,c(upgn_cir_gnset$term_size[i],
                                              upgn_cir_gnset$p_value[i],upgn_cir_gnset$term_id[i]
                                              ,upgn_cir_gnset$term_name[i],upgn_cir_gnset$query_size[i]
                                              ,upgn_cir_gnset$intersection_size[i],upgn_cir_gnset$intersection[i]))
  }
}

colnames(upgn_cir_gnset2 )=c('term_size',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
upgn_cir_gnset2=upgn_cir_gnset2[order(as.numeric(upgn_cir_gnset2$p_value)),]

write.csv(upgn_cir_gnset2, file = "upgn_cir_gnset_reanno.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)



#這個元資料的 GSEA跑完ontol lip 後 在作一個top 100 node 的 GSEA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

upgn_mrna_gnset_100top=gprofiler2::gost(freq_mrna_up_100top$Var1, correction_method="fdr",
                                 sources =c("GO:BP","GO:MF"),evcodes = TRUE)
dwgn_mrna_gnset_100top=gprofiler2::gost(freq_mrna_dw_100top$Var1, correction_method="fdr",
                                 sources =c("GO:BP","GO:MF"),evcodes = TRUE)

upgn_cir_gnset_100top=gprofiler2::gost(freq_cir_up_100top$Var1, correction_method="fdr",
                                 sources =c("GO:BP","GO:MF"),evcodes = TRUE)
dwgn_cir_gnset_100top=gprofiler2::gost(freq_cir_dw_100top$Var1, correction_method="fdr",
                                sources =c("GO:BP","GO:MF"),evcodes = TRUE)

upgn_mrna_gnset_100top=upgn_mrna_gnset_100top$result
dwgn_mrna_gnset_100top=dwgn_mrna_gnset_100top$result
upgn_cir_gnset_100top=upgn_cir_gnset_100top$result
dwgn_cir_gnset_100top=dwgn_cir_gnset_100top$result


upgn_mrna_gnset_100top2<- data.frame()
for (i in 1:length(upgn_mrna_gnset_100top$term_size)){
  if(upgn_mrna_gnset_100top$term_size[i]<1000 && upgn_mrna_gnset_100top$term_size[i]>10){
    print(i)
    print(upgn_mrna_gnset_100top$term_size[i])
    upgn_mrna_gnset_100top2=rbind(upgn_mrna_gnset_100top2,c(upgn_mrna_gnset_100top$term_size[i],
                                            upgn_mrna_gnset_100top$p_value[i],upgn_mrna_gnset_100top$term_id[i]
                                            ,upgn_mrna_gnset_100top$term_name[i],upgn_mrna_gnset_100top$query_size[i]
                                            ,upgn_mrna_gnset_100top$intersection_size[i],upgn_mrna_gnset_100top$intersection[i]))
  }
}

colnames(upgn_mrna_gnset_100top2 )=c('term_size',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
upgn_mrna_gnset_100top2=upgn_mrna_gnset_100top2[order(as.numeric(upgn_mrna_gnset_100top2$p_value)),]

write.csv(upgn_mrna_gnset_100top2, file = "upgn_mrna_gnset_100top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


dwgn_mrna_gnset_100top2<- data.frame()
for (i in 1:length(dwgn_mrna_gnset_100top$term_size)){
  if(dwgn_mrna_gnset_100top$term_size[i]<1000 && dwgn_mrna_gnset_100top$term_size[i]>10){
    print(i)
    print(dwgn_mrna_gnset_100top$term_size[i])
    dwgn_mrna_gnset_100top2=rbind(dwgn_mrna_gnset_100top2,c(dwgn_mrna_gnset_100top$term_size[i],
                                                            dwgn_mrna_gnset_100top$p_value[i],dwgn_mrna_gnset_100top$term_id[i]
                                                            ,dwgn_mrna_gnset_100top$term_name[i],dwgn_mrna_gnset_100top$query_size[i]
                                                            ,dwgn_mrna_gnset_100top$intersection_size[i],dwgn_mrna_gnset_100top$intersection[i]))
  }
}

colnames(dwgn_mrna_gnset_100top2 )=c('term_size',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
dwgn_mrna_gnset_100top2=dwgn_mrna_gnset_100top2[order(as.numeric(dwgn_mrna_gnset_100top2$p_value)),]

write.csv(dwgn_mrna_gnset_100top2, file = "dwgn_mrna_gnset_100top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#cir
upgn_cir_gnset_100top2<- data.frame()
for (i in 1:length(upgn_cir_gnset_100top$term_size)){
  if(upgn_cir_gnset_100top$term_size[i]<1000 && upgn_cir_gnset_100top$term_size[i]>10){
    print(i)
    print(upgn_cir_gnset_100top$term_size[i])
    upgn_cir_gnset_100top2=rbind(upgn_cir_gnset_100top2,c(upgn_cir_gnset_100top$term_size[i],
                                                            upgn_cir_gnset_100top$p_value[i],upgn_cir_gnset_100top$term_id[i]
                                                            ,upgn_cir_gnset_100top$term_name[i],upgn_cir_gnset_100top$query_size[i]
                                                            ,upgn_cir_gnset_100top$intersection_size[i],upgn_cir_gnset_100top$intersection[i]))
  }
}

colnames(upgn_cir_gnset_100top2 )=c('term_size',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
upgn_cir_gnset_100top2=upgn_cir_gnset_100top2[order(as.numeric(upgn_cir_gnset_100top2$p_value)),]

write.csv(upgn_cir_gnset_100top2, file = "upgn_cir_gnset_100top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


dwgn_cir_gnset_100top2<- data.frame()
for (i in 1:length(dwgn_cir_gnset_100top$term_size)){
  if(dwgn_cir_gnset_100top$term_size[i]<1000 && dwgn_cir_gnset_100top$term_size[i]>10){
    print(i)
    print(dwgn_cir_gnset_100top$term_size[i])
    dwgn_cir_gnset_100top2=rbind(dwgn_cir_gnset_100top2,c(dwgn_cir_gnset_100top$term_size[i],
                                                          dwgn_cir_gnset_100top$p_value[i],dwgn_cir_gnset_100top$term_id[i]
                                                          ,dwgn_cir_gnset_100top$term_name[i],dwgn_cir_gnset_100top$query_size[i]
                                                          ,dwgn_cir_gnset_100top$intersection_size[i],dwgn_cir_gnset_100top$intersection[i]))
  }
}

colnames(dwgn_cir_gnset_100top2 )=c('term_size',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
dwgn_cir_gnset_100top2=dwgn_cir_gnset_100top2[order(as.numeric(dwgn_cir_gnset_100top2$p_value)),]

write.csv(dwgn_cir_gnset_100top2, file = "dwgn_cir_gnset_100top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


#看原本900組的 rna 根 top 100 node 的共同gnset有幾個
# make geneset up-up list and dw-dw list
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
upgn_mrna_gnset_100top2 = read.table(file = "upgn_mrna_gnset_100top.csv", header = TRUE, sep = ",")
dwgn_mrna_gnset_100top2 = read.table(file = "dwgn_mrna_gnset_100top.csv", header = TRUE, sep = ",")
upgn_cir_gnset_100top2 = read.table(file = "upgn_cir_gnset_100top.csv", header = TRUE, sep = ",")
dwgn_cir_gnset_100top2 = read.table(file = "dwgn_cir_gnset_100top.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/ontol_up")
upgn_mrna_gnset2= read.table(file = "upgn_mrna_gnset_reanno.csv", header = TRUE, sep = ",")
upgn_cir_gnset2= read.table(file = "upgn_cir_gnset_reanno.csv", header = TRUE, sep = ",")
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/ontol_dw")
dwgn_mrna_gnset2= read.table(file = "dwgn_mrna_gnset_reanno.csv", header = TRUE, sep = ",")
dwgn_cir_gnset2= read.table(file = "dwgn_cir_gnset(ori).csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")


up_up_mrna= list("up_mRNA_reanno"=upgn_mrna_gnset2$term_ID,"up_mRNA_top100"=upgn_mrna_gnset_100top2$term_ID)
up_up_cir= list("up_cir_reanno"=upgn_cir_gnset2$term_ID,"up_cir_top100"=upgn_cir_gnset_100top2$term_ID)

dw_dw_mrna= list("dw_mRNA_reanno"=dwgn_mrna_gnset2$term_ID,"dw_mRNA_top100"=dwgn_mrna_gnset_100top2$term_ID)
dw_dw_cir= list("dw_cir_reanno"=dwgn_cir_gnset2$term_ID,"dw_cir_top100"=dwgn_cir_gnset_100top2$term_ID)


# venndiagram
venn.diagram(up_up_mrna, filename = 'up_up_mrna_gnset.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'),col = c('red', 'blue'), cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif',margin=0.05,cat.pos=0)
venn.diagram(up_up_cir, filename = 'up_up_cir_gnset.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'),col = c('red', 'blue'), cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif',margin=0.05,cat.pos=0)
venn.diagram(dw_dw_mrna, filename = 'dw_dw_mrna_gnset.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'),col = c('red', 'blue'), cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif',margin=0.05,cat.pos=0) #,cat.dist=0.03
venn.diagram(dw_dw_cir, filename = 'dw_dw_cir_gnset.png', imagetype = 'png', fill = c('red', 'blue'), alpha = 0.50, 
             cat.col = c('red', 'blue'),col = c('red', 'blue'), cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif',margin=0.05,cat.pos=0)

# 取得 common gnset 1102 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

up_up_mrna_gnset= merge(x=upgn_mrna_gnset2, y=upgn_mrna_gnset_100top2, by.x=c("term_ID"), by.y = c("term_ID"))
up_up_cir_gnset= merge(x=upgn_cir_gnset2, y=upgn_cir_gnset_100top2, by.x=c("term_ID"), by.y = c("term_ID"))
dw_dw_mrna_gnset=  merge(x=dwgn_mrna_gnset2, y=dwgn_mrna_gnset_100top2, by.x=c("term_ID"), by.y = c("term_ID"))
dw_dw_cir_gnset= merge(x=dwgn_cir_gnset2, y=dwgn_cir_gnset_100top2, by.x=c("term_ID"), by.y = c("term_ID"))

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")

annodata= read.csv(file = "AllAnnoData2.csv", header = TRUE, sep = "\t",quote = "", stringsAsFactors=FALSE )

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/common")
up_up_mrna_gnset2= merge(x=up_up_mrna_gnset, y=annodata, by.x=c("term_name.x"), by.y = c("term_name"))
up_up_cir_gnset2= merge(x=up_up_cir_gnset, y=annodata, by.x=c("term_name.x"), by.y = c("term_name"))
dw_dw_mrna_gnset2= merge(x=dw_dw_mrna_gnset, y=annodata, by.x=c("term_name.x"), by.y = c("term_name"))
dw_dw_cir_gnset2= merge(x=dw_dw_cir_gnset, y=annodata, by.x=c("term_name.x"), by.y = c("term_name"))

#remove duplicate
dw_dw_mrna_gnset2= dw_dw_mrna_gnset2 %>% distinct()

#在首動 作表格  
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/common")

write.csv(up_up_mrna_gnset2, file = "up_up_mrna_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(up_up_cir_gnset2, file = "up_up_cir_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dw_dw_mrna_gnset2, file = "dw_dw_mrna_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dw_dw_cir_gnset2, file = "dw_dw_cir_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# 加上顏色 作長條圖
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")
color_code= read.csv(file = "color_code.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/common")
dw_dw_mrna_gnset= read.csv(file = "dw_dw_mrna_gnset22.csv", header = TRUE, sep = ",")

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

# 取得 common gnset 1102 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# 10/20 2023    group reanno 的到8種種類
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")

annodata= read.csv(file = "AllAnnoData.csv", header = TRUE, sep = ",",quote = "", stringsAsFactors=FALSE )
length(intersect(upgn_mrna_gnset2$term_name,annodata$term_name))
length(intersect(dwgn_mrna_gnset2$term_name,annodata$term_name))
length(intersect(upgn_cir_gnset2$term_name,annodata$term_name))
length(intersect(dwgn_cir_gnset2$term_name,annodata$term_name))
upgn_mrna_gnset3= upgn_mrna_gnset2[upgn_mrna_gnset2$term_name %in% annodata$term_name,]
upgn_mrna_gnset2_not= upgn_mrna_gnset2[!upgn_mrna_gnset2$term_name %in% annodata$term_name,]
dwgn_mrna_gnset2_not= dwgn_mrna_gnset2[!dwgn_mrna_gnset2$term_name %in% annodata$term_name,]
upgn_cir_gnset2_not= upgn_cir_gnset2[!upgn_cir_gnset2$term_name %in% annodata$term_name,]
dwgn_cir_gnset2_not= dwgn_cir_gnset2[!dwgn_cir_gnset2$term_name %in% annodata$term_name,]

write.csv(upgn_mrna_gnset2_not, file = "upgn_mrna_gnset2_not.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_mrna_gnset2_not, file = "dwgn_mrna_gnset2_not.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(upgn_cir_gnset2_not, file = "upgn_cir_gnset2_not.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_cir_gnset2_not, file = "dwgn_cir_gnset2_not.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# 這個讀ann檔案 怪怪的 用python作  mrna序列且mi-mrna結果.ipynb

#存到 /home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reanno_group
# upgn_mrna_gnset_reanno_group.csv

# 在來作pie chart 
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reanno_group")
upgn_mrna_gnset3= read.csv(file = "upgn_mrna_gnset_reanno_group.csv", header = TRUE, sep = ",")
upgn_cir_gnset3= read.csv(file = "upgn_cir_gnset_reanno_group.csv", header = TRUE, sep = ",")
dwgn_mrna_gnset3= read.csv(file = "dwgn_mrna_gnset_reanno_group.csv", header = TRUE, sep = ",")
dwgn_cir_gnset3= read.csv(file = "dwgn_cir_gnset(ori)_group.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy")
color_code= read.csv(file = "color_code.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reanno_group")

# dwgn_mrna_gnset3 又部能讀檔案 去用python作   mrna序列且mi-mrna結果.ipynb
# 發現用read csv就可以讀了
upgn_mrna_gnset_group= upgn_mrna_gnset3 %>%
  group_by(Group) %>%
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

upgn_mrna_gnset_group2= merge(x=upgn_mrna_gnset_group, y=color_code, by.x=c("Group"), by.y = c("group"))

# pie chart  
# geom_bar(stat = "identity") coord_polar(theta = "y", start=0) 這是一定要的

ggplot(upgn_mrna_gnset_group, aes(x="" ,y=perc, fill=Group)) +
  geom_bar(stat = "identity")+
  coord_polar(theta = "y", start=0)+
  geom_col(color = "white") +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  theme_void()
  
# 加上所用的顏色
ggplot(upgn_mrna_gnset_group, aes(x="" ,y=perc, fill=Group)) +
  geom_bar(stat = "identity")+
  coord_polar(theta = "y", start=0)+
  geom_col(color = "white") +
  scale_fill_manual(values = upgn_mrna_gnset_group2$code) +
  labs(x = "", y = "", title = "upgn_mrna_gnset_group \n",
       fill = "Group") + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_text(hjust = 0.5, face="bold", size = 10))+
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  theme_void()


# 其他組
upgn_cir_gnset_group= upgn_cir_gnset3 %>%
  group_by(Group) %>%
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

upgn_cir_gnset_group2= merge(x=upgn_cir_gnset_group, y=color_code, by.x=c("Group"), by.y = c("group"))

dwgn_mrna_gnset_group= dwgn_mrna_gnset3 %>%
  group_by(Group) %>%
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

dwgn_mrna_gnset_group2= merge(x=dwgn_mrna_gnset_group, y=color_code, by.x=c("Group"), by.y = c("group"))

dwgn_cir_gnset_group= dwgn_cir_gnset3 %>%
  group_by(Group) %>%
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

dwgn_cir_gnset_group2= merge(x=dwgn_cir_gnset_group, y=color_code, by.x=c("Group"), by.y = c("group"))

#畫圖
ggplot(upgn_cir_gnset_group, aes(x="" ,y=perc, fill=Group)) +
  geom_bar(stat = "identity")+
  coord_polar(theta = "y", start=0)+
  geom_col(color = "white") +
  scale_fill_manual(values = upgn_cir_gnset_group2$code) +
  labs(x = "", y = "", title = "upgn_cir_gnset_group",
       fill = "Group") + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_text(hjust = 0.5, face="bold", size = 10))+
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  theme_void()

ggplot(dwgn_mrna_gnset_group, aes(x="" ,y=perc, fill=Group)) +
  geom_bar(stat = "identity")+
  coord_polar(theta = "y", start=0)+
  geom_col(color = "white") +
  scale_fill_manual(values = dwgn_mrna_gnset_group2$code) +
  labs(x = "", y = "", title = "dwgn_mrna_gnset_group",
       fill = "Group") + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_text(hjust = 0.5, face="bold", size = 10))+
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  theme_void()

ggplot(dwgn_cir_gnset_group, aes(x="" ,y=perc, fill=Group)) +
  geom_bar(stat = "identity")+
  coord_polar(theta = "y", start=0)+
  geom_col(color = "white") +
  scale_fill_manual(values = dwgn_cir_gnset_group2$code) +
  labs(x = "", y = "", title = "dwgn_cir_gnset_group",
       fill = "Group") + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_text(hjust = 0.5, face="bold", size = 10))+
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  theme_void()



#作 bar plot 看看
#  upgn_mrna_gnset3  g1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
category_counts <- table(upgn_mrna_gnset3$Group)
category_counts <- sort(category_counts, decreasing = TRUE)
category_counts <- data.frame(category_counts)
s <- sum(category_counts$Freq)
category_counts$Freq <- category_counts$Freq / sum(category_counts$Freq)*100

category_counts2= merge(x=category_counts, y=color_code, by.x=c("Var1"), by.y = c("group"))
category_counts2=category_counts2[order(category_counts2$Freq, decreasing = TRUE),]

ggplot(category_counts2, aes(x = Var1, y = Freq, fill=Var1)) +
  geom_bar(stat = "identity") +
  labs(title = paste("upgn_mrna_gnset \n","n =",as.character(s)), x = "Group", y = "percentage") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1,size=11))+
  scale_fill_manual(values = category_counts2$code)+
  theme(legend.position="none")

g1= ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=round(as.numeric(Freq), digits = 1)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip(ylim=c(0,65))  + #把圖轉至
  labs( x = "", y = "",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic(base_size = 15)  + #把所有原本字體同時設定
  theme()+
  theme(legend.position = "none") #要有這個把legend刪掉



#  upgn_cir_gnset3  g2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>
category_counts <- table(upgn_cir_gnset3$Group)
category_counts <- sort(category_counts, decreasing = TRUE)
category_counts <- data.frame(category_counts)
s <- sum(category_counts$Freq)
category_counts$Freq <- category_counts$Freq / sum(category_counts$Freq)*100

category_counts2= merge(x=category_counts, y=color_code, by.x=c("Var1"), by.y = c("group"))
category_counts2=category_counts2[order(category_counts2$Freq, decreasing = TRUE),]

ggplot(category_counts2, aes(x = Var1, y = Freq, fill=Var1)) +
  geom_bar(stat = "identity") +
  labs(title = paste("upgn_cir_gnset \n","n =",as.character(s)), x = "Group", y = "percentage") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1,size=11))

# format + round 可以把60.0 秀出來
g2= ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=format(round(as.numeric(Freq),digits = 1), nsmall = 1)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip(ylim=c(0,65))  + #把圖轉至
  labs( x = "", y = "",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic(base_size = 15)  + #把所有原本字體同時設定
  theme()+
  theme(legend.position = "none") #要有這個把legend刪掉



#dwgn_mrna_gnset3  g3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#這邊嘗試把圖做成恆的 加上數字在塗上  還有想要剪圖 把scale弄小一點 
category_counts <- table(dwgn_mrna_gnset3$Group)
category_counts <- sort(category_counts, decreasing = TRUE)
category_counts <- data.frame(category_counts[1:12])
s <- sum(category_counts$Freq)
category_counts$Freq <- category_counts$Freq / sum(category_counts$Freq)*100


category_counts2= merge(x=category_counts, y=color_code, by.x=c("Var1"), by.y = c("group"))
category_counts2=category_counts2[order(category_counts2$Freq, decreasing = TRUE),]

#做好圖 加上數字
ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=round(as.numeric(Freq), digits = 3)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip()  + #把圖轉至
  labs(title = paste("Downregulated mRNA geneset  \n","n =",as.character(s)), x = "", y = "percentage",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic(base_size = 13)  + #把所有原本字體同時設定
  theme(axis.text.x = element_text(hjust = 1,size=11),
        axis.text.y = element_text(hjust = 1,size=11))+
  theme(legend.position = "none") #要有這個把legend刪掉

#存g3 就能家在一起
g3= ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=round(as.numeric(Freq), digits = 1)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip(ylim=c(0,65))  + #把圖轉至
  labs( x = "", y = "",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic(base_size = 15)  + #把所有原本字體同時設定
  theme()+
  theme(legend.position = "none") #要有這個把legend刪掉


#現在要分途 把圖剪裁 縮小sacle >>>>>>>>>>>>>>>>>>>>>>
#需要 gridextra
library(gridExtra)

g1= ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=round(as.numeric(Freq), digits = 1)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip(ylim=c(40,55))+
  labs( x = "", y = "percentage",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic(base_size = 13)  + #把所有原本字體同時設定
  theme(axis.text.y = element_blank(),
        #axis.text.x = element_blank(), 這個會有x軸的標籤 scale
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) + 
  theme(legend.position = "none") #要有這個把legend刪掉


g2=  ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=round(as.numeric(Freq), digits = 1)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip(ylim=c(0,15))+
  labs( x = "", y = "percentage",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  
  theme_classic(base_size = 13)  + #把所有原本字體同時設定
  theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  theme(legend.position = "none") #要有這個把legend刪掉


grid.arrange(g2,g1, heights=c(3/4, 3/4),
             widths=c(3/4, 1/4), ncol=2, nrow=2)

#雖然做到 簡短圖 但是不夠好 就算了
#現在要分途 把圖剪裁 縮小sacle >>>>>>>>>>>>>>>>>>>>>>

  

#  dwgn_cir_gnset3  g4 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
category_counts <- table(dwgn_cir_gnset3$Group)
category_counts <- sort(category_counts, decreasing = TRUE)
category_counts <- data.frame(category_counts)
s <- sum(category_counts$Freq)
category_counts$Freq <- category_counts$Freq / sum(category_counts$Freq)*100 

category_counts2= merge(x=category_counts, y=color_code, by.x=c("Var1"), by.y = c("group"))
category_counts2=category_counts2[order(category_counts2$Freq, decreasing = TRUE),]

ggplot(category_counts2, aes(x = Var1, y = Freq, fill=Var1)) +
  geom_bar(stat = "identity") +
  labs(title = paste("dwgn_cir_gnset \n","n =",as.character(s)), x = "Group", y = "percentage") +
    scale_fill_manual(values = category_counts2$code)+
  theme(legend.position="none")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1,size=16))

g4= ggplot(category_counts2, aes(x = reorder(Var1, Freq), y = Freq, fill=Var1), show.legend = FALSE) + #要有這個把legend刪掉
  geom_text(aes(label=format(round(as.numeric(Freq), digits = 1),nsmall=1)), hjust=-0.5, colour = "black", position = "dodge") +
  geom_bar(aes(fill=Var1), stat = "identity", position = position_dodge(width = .25))+
  coord_flip(ylim=c(0,65))  + #把圖轉至
  labs( x = "", y = "",fill = "") +
  scale_fill_manual(values = category_counts2$code)+
  theme_classic(base_size = 15)  + #把所有原本字體同時設定
  theme()+
  theme(legend.position = "none") #要有這個把legend刪掉


#把所有圖組合 >>>>>>>>>>>>>>>>>>>>>>>>>>>
# cowplot import
plot_grid(g1,g3,g2,g4,labels = c('A', 'B', 'C', 'D'), label_size = 19,align="hv",ncol=2, nrow=2)

#把所有圖組合 >>>>>>>>>>>>>>>>>>>>>>>>>>>


#作長條圖 前10 gnset  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
upgn_mrna_gnset3_top10= upgn_mrna_gnset3[1:10,] 
upgn_mrna_gnset3_top10$log_p_value= -log10(upgn_mrna_gnset3_top10$p_value)


ggplot(upgn_mrna_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#4c6cb3")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 25))+
  labs(title = "", x = "", y = "p value (-log10)")

g1= ggplot(upgn_mrna_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#4c6cb3")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 18))+
  labs(title = "", x = "", y = "p value (-log10)")

# up cir gnset top 10
upgn_cir_gnset3_top10= upgn_cir_gnset3[1:10,] 
upgn_cir_gnset3_top10$log_p_value= -log10(upgn_cir_gnset3_top10$p_value)

ggplot(upgn_cir_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#c97586")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 25))+
  labs(title = "", x = "", y = "p value (-log10)")

g2= ggplot(upgn_cir_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#c97586")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 18))+
  labs(title = "", x = "", y = "p value (-log10)")

# dw mrna  top 10
dwgn_mrna_gnset3_top10= dwgn_mrna_gnset3[1:10,] 
dwgn_mrna_gnset3_top10$log_p_value= -log10(dwgn_mrna_gnset3_top10$p_value)


ggplot(dwgn_mrna_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#4c6cb3")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 25))+
  labs(title = "", x = "", y = "p value (-log10)")

g3= ggplot(dwgn_mrna_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#4c6cb3")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 18))+
  labs(title = "", x = "", y = "p value (-log10)")


# dw cir  gnset top 10
dwgn_cir_gnset3_top10= dwgn_cir_gnset3[1:10,] 
dwgn_cir_gnset3_top10$log_p_value= -log10(dwgn_cir_gnset3_top10$p_value)

ggplot(dwgn_cir_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#c97586")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 25))+
  labs(title = "", x = "", y = "p value (-log10)")

g4= ggplot(dwgn_cir_gnset3_top10, aes(x =reorder(term_name, log_p_value) , y = log_p_value)) +
  geom_bar(stat = "identity", fill=c("#c97586")) +
  coord_flip() +
  theme(legend.position="none")+
  scale_y_log10()+
  theme_classic()+
  theme(text = element_text(size = 18))+
  labs(title = "", x = "", y = "p value (-log10)")


#把所有圖組合 >>>>>>>>>>>>>>>>>>>>>>>>>>>
# cowplot import
library(cowplot)
plot_grid(g1,g3,g2,g4,labels = c('A', 'B', 'C', 'D'),label_size = 20, align="hv",ncol=2, nrow=2)

#把所有圖組合 >>>>>>>>>>>>>>>>>>>>>>>>>>>


# 11 02  把mrna的 gsea 的metabolic gnset  cell cycle取出來一個table

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA/genesetRedundancy/reanno_group")
upgn_mrna_gnset3= read.csv(file = "upgn_mrna_gnset_reanno_group.csv", header = TRUE, sep = ",")
dwgn_mrna_gnset3= read.csv(file = "dwgn_mrna_gnset_reanno_group.csv", header = TRUE, sep = ",")


upgn_mrna_gnset_group_sub= upgn_mrna_gnset3 %>%
  group_by(Group) %>%
  filter(Group == "metabolic process") 

dwgn_mrna_gnset_group_sub= dwgn_mrna_gnset3 %>%
  group_by(Group) %>%
  filter(Group == "cell cycle") 

write.csv(upgn_mrna_gnset_group_sub, file = "upgn_mrna_gnset_metabolic.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_mrna_gnset_group_sub, file = "dwgn_mrna_gnset_cellcycle.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

