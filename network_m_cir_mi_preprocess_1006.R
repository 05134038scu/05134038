# pearson cor of mrna cir
library(doParallel)
library(foreach)


setwd("/home/yuanpeng/yuan/R sarsCoV2/result/redo_0920")
mrna_fr<- read.table(file = "updwgn_fr_mrna3.csv", header = TRUE, sep = ",",row.names=1)

setwd("/home/yuanpeng/yuan/R sarsCoV2/cirRNA_ensg/new2/redo_0920")
cir_fr<- read.table(file = "updwgn_fr_cir3.csv", header = TRUE, sep = ",",row.names=1)

setwd("/home/yuanpeng/yuan/R sarsCoV2/mirna_sequen")
mir_fr<- read.table(file = "mirna_fr3.csv", header = TRUE, sep = ",",row.names=1)

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/pearson")

mrna_fr =mrna_fr[,-c(12:14)]
cir_fr=cir_fr[,-c(12:14)]

colnames(mrna_fr)=c("0_r2_notinf","0_r3_notinf","3_r3_inf","6_r1_inf","6_r2_inf","6_r3_inf",
                    "12_r1_inf","12_r3_inf","24_r1_inf","24_r2_inf","24_r3_inf" )
colnames(cir_fr)=c("0_r2_notinf","0_r3_notinf","3_r3_inf","6_r1_inf","6_r2_inf","6_r3_inf",
                    "12_r1_inf","12_r3_inf","24_r1_inf","24_r2_inf","24_r3_inf")

colnames(mir_fr)=c("0_r2_notinf","0_r3_notinf","3_r3_inf","6_r2_inf","6_r3_inf",
                   "12_r3_inf","24_r2_inf","24_r3_inf")

mrna_fr2=t(mrna_fr)
cir_fr2=t(cir_fr)
mir_fr2=t(mir_fr)
# log data 才能作相關性
mrna_fr3=log2(mrna_fr2+1)
cir_fr3=log2(cir_fr2+1)
mir_fr3=log2(mir_fr2+1)


#origin code
mrna_cir_pear.cor_all_res=data.frame()
mrna_cir_pear.cor_all=data.frame()

for (i in  1:ncol(mrna_fr2)) {
  for(j in 1:ncol(cir_fr2)){
    mrna_cir_pear.cor_all_res = cor.test(mrna_fr2[,i], cir_fr2[,j], alternative="two.side", method="pearson", conf.level=0.95)
    mrna_cir_pear.cor_all=rbind(mrna_cir_pear.cor_all,c(mrna_cir_pear.cor_all_res$p.value,mrna_cir_pear.cor_all_res$estimate,colnames(mrna_fr2)[i],colnames(cir_fr2)[j]))
  }
}

# use parallel 跑 對兩組資料互相跑 test >>>>>>>>>>>>>>>>>>>>
aa=foreach(i=1:100, .combine="rbind") %do% rnorm(4)
x=foreach(i=1:100, .combine="rbind") %do% rnorm(6)
num_cores <- 4  # Change this to the number of CPU cores you want to use

cl <- makeCluster(num_cores)
registerDoParallel(cl)
results_list <- foreach(i = 1:ncol(aa), .combine = rbind) %:%
  foreach(j = 1:ncol(x), .combine = rbind) %dopar% {
    res <- cor.test(aa[, i], x[, j], alternative = "two.side", method = "pearson", conf.level = 0.95)
    # Add results to the data frame
    data.frame(
      p.value = res$p.value,
      estimate = res$estimate
      
    )
  }
stopCluster(cl)

# 所以時用上
results_list=c()
num_cores <- 16
cl <- makeCluster(num_cores)
registerDoParallel(cl)
results_list <- foreach(i = 1:ncol(mrna_fr3), .combine = rbind) %:%
  foreach(j = 1:ncol(cir_fr3), .combine = rbind) %dopar% {
    res <- cor.test(mrna_fr3[, i], cir_fr3[, j], alternative = "two.side", method = "pearson", conf.level = 0.95)
    # Add results to the data frame
    data.frame(
      p.value = res$p.value,estimate = res$estimate,colnames(mrna_fr3)[i],colnames(cir_fr3)[j]
    )
  }
stopCluster(cl)


colnames(results_list)= c("pvalue","pearson correlation", "mRNA_ENSG", "cir_ENST")
results_list2=results_list[results_list$pvalue<0.05 & results_list$`pearson correlation`> 0.7,]

write.csv(results_list2, file = "mrna_cir_pear.cor_0.05p.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)


# use parallel 重複跑100次 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cl <- makeCluster(num_cores)
registerDoParallel(cl)
num_iterations <- 100
#emp <- x[sample(nrow(x)), ]
results_list <- foreach(iteration = 1:num_iterations, .combine = rbind) %dopar% {
  pear.cor <- data.frame()
  temp <- x[sample(nrow(x)), ]
  for (i in 1:ncol(aa)) {
    for(j in 1:ncol(x)){
      res <- cor.test(aa[, i], x[, j], alternative = "two.side", method = "pearson", conf.level = 0.95)
      pear.cor <- rbind(pear.cor, c(res$p.value, res$estimate,aa[,]))
    }
  }
  
  return(pear.cor)  # Return the data frame for this iteration
}
stopCluster(cl)

# use parallel 重複跑100次 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# mrna mirna perason corr, cir mi perason corr
# remove row of mrna cir
mrna_fr3=mrna_fr3[-c(4,7,9),]
cir_fr3=cir_fr3[-c(4,7,9),]

# mrna mirna perason corr
results_list=c()
num_cores <- 16
cl <- makeCluster(num_cores)
registerDoParallel(cl)
results_list <- foreach(i = 1:ncol(mrna_fr3), .combine = rbind) %:%
  foreach(j = 1:ncol(mir_fr3), .combine = rbind) %dopar% {
    res <- cor.test(mrna_fr3[, i], mir_fr3[, j], alternative = "two.side", method = "pearson", conf.level = 0.95)
    # Add results to the data frame
    data.frame(
      p.value = res$p.value,estimate = res$estimate,colnames(mrna_fr3)[i],colnames(mir_fr3)[j]
    )
  }
stopCluster(cl)


colnames(results_list)= c("pvalue","pearson correlation", "mRNA_ENSG", "mir_ENST")
results_list2=results_list[results_list$pvalue<0.05 & results_list$`pearson correlation`< -0.7,]

write.csv(results_list2, file = "mrna_mi_pear.cor_0.05p.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)

# cir mi pearson corr

results_list=c()
num_cores <- 16
cl <- makeCluster(num_cores)
registerDoParallel(cl)
results_list <- foreach(i = 1:ncol(cir_fr3), .combine = rbind) %:%
  foreach(j = 1:ncol(mir_fr3), .combine = rbind) %dopar% {
    res <- cor.test(cir_fr3[, i], mir_fr3[, j], alternative = "two.side", method = "pearson", conf.level = 0.95)
    # Add results to the data frame
    data.frame(
      p.value = res$p.value,estimate = res$estimate,colnames(cir_fr3)[i],colnames(mir_fr3)[j]
    )
  }
stopCluster(cl)


colnames(results_list)= c("pvalue","pearson correlation", "cir_ENST", "mir_ENST")
results_list2=results_list[results_list$pvalue<0.05 & results_list$`pearson correlation`< -0.7,]

write.csv(results_list2, file = "cir_mi_pear.cor_0.05p.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)

# get unique mrna cir mi ENSG 
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/pearson")

mrna_cir_pear<- read.table(file = "mrna_cir_pear.cor_0.05p.csv", header = TRUE, sep = ",")
mrna_mi_pear<- read.table(file = "mrna_mi_pear.cor_0.05p.csv", header = TRUE, sep = ",")
cir_mi_pear<- read.table(file = "cir_mi_pear.cor_0.05p.csv", header = TRUE, sep = ",")

mrna_com= intersect(mrna_cir_pear$mRNA_ENSG,mrna_mi_pear$mRNA_ENSG)
cir_com= intersect(mrna_cir_pear$cir_ENST,cir_mi_pear$cir_ENST)
mi_com= intersect(cir_mi_pear$mir_ENST,mrna_mi_pear$mir_ENST)




# 轉 mrna gn id
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

mrna_ID=biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name'),
                      filters = 'ensembl_gene_id', 
                      values = mrna_com,
                      mart = mart)

# 轉 cir gn id
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

cir_ID=biomaRt::getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id','hgnc_symbol','external_gene_name'),
                      filters = 'ensembl_transcript_id', 
                      values = cir_com,
                      mart = mart)


mrna_fr_com=mrna_fr[row.names(mrna_fr) %in% mrna_com, ]
cir_fr_com=cir_fr[row.names(cir_fr) %in% cir_com, ]


setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/seq")
mi_pair<- read.table(file = "mi_pair.txt", header = TRUE, sep = "\t")
colnames(mi_pair)=c('mi_enst','mi_id')

mi_pair2 <- mi_pair[mi_pair$mi_enst %in% mi_com , ]
length(intersect(mi_com, mi_pair$mi_enst))

mi_fr_com=mir_fr[row.names(mir_fr) %in% mi_pair2$mi_enst, ]

#存 剩下的共同 mrna and cir mi fr
write.csv(mrna_fr_com, file = "mrna_com_cir_mi_fr.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
write.csv(cir_fr_com, file = "cir_com_mrna_mi_fr.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
write.csv(mi_fr_com, file = "mir_com_mrna_cir_fr.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)


write.csv(mi_pair2, file = "mi_pair2.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)

#把這些有共同的 mrna cir mi 取得名字去根miranda結果比對
# mi 已經作了
setwd("/home/yuanpeng/yuan/R sarsCoV2/result/redo_0920")
mrna_fr<- read.table(file = "updwgn_fr_mrna3.csv", header = TRUE, sep = ",",row.names=1)

setwd("/home/yuanpeng/yuan/R sarsCoV2/cirRNA_ensg/new2/redo_0920")
cir_fr<- read.table(file = "updwgn_fr_cir3.csv", header = TRUE, sep = ",",row.names=1)
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/pearson")


colnames(mrna_fr)=c("0_r2_notinf","0_r3_notinf","3_r3_inf","6_r1_inf","6_r2_inf","6_r3_inf",
                    "12_r1_inf","12_r3_inf","24_r1_inf","24_r2_inf","24_r3_inf" ,"ENSG","hgnc","gn_name")
colnames(cir_fr)=c("0_r2_notinf","0_r3_notinf","3_r3_inf","6_r1_inf","6_r2_inf","6_r3_inf",
                   "12_r1_inf","12_r3_inf","24_r1_inf","24_r2_inf","24_r3_inf","ENST","ENSG","gn_name")


mrna_fr_com=mrna_fr[row.names(mrna_fr) %in% mrna_com, ]
cir_fr_com=cir_fr[row.names(cir_fr) %in% cir_com, ]

write.csv(mrna_fr_com, file = "mrna_com_cir_mi_fr2.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)
write.csv(cir_fr_com, file = "cir_com_mrna_mi_fr2.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)



#>>>>>>>>>>>>>>>>>>>10/01  new network 前至步驟 作 network table 根att list>>>>>>>>>>>>>>>>>>>>>>
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2")

mrna_fin=read.table(file = "mrna_mi_cir_predict_fin.csv", header = TRUE, sep = ",")
cir_fin=read.table(file = "cir_mi_cir_predict_fin.csv", header = TRUE, sep = ",")

#first group sub gn
cir_fin_sub_gn <- cir_fin %>%
  group_by(miRNA, gene) %>%
  count() %>%
  rename(node_degree = n)


mrna_fin_sub_gn <- mrna_fin %>%
  group_by(miRNA, mRNA) %>%
  count() %>%
  rename(node_degree = n)

write.csv(cir_fin_sub_gn, file = "cir_fin_sub_gn.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(mrna_fin_sub_gn, file = "mrna_fin_sub_gn.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# group max score 

cir_fin_max_score <- cir_fin %>%
  group_by(miRNA, gene) %>%
  arrange(desc(score)) %>%
  slice(1)


mrna_fin_max_score <- mrna_fin %>%
  group_by(miRNA, mRNA) %>%
  arrange(desc(score)) %>%
  slice(1)

write.csv(cir_fin_max_score, file = "cir_fin_max_score.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(mrna_fin_max_score , file = "mrna_fin_max_score.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#merge 2 table

cir_fin_merged <- merge(cir_fin_sub_gn, cir_fin_max_score, by = c("miRNA", "gene"), all.x = TRUE)
mrna_fin_merged <- merge(mrna_fin_sub_gn, mrna_fin_max_score, by = c("miRNA", "mRNA"), all.x = TRUE)

write.csv(cir_fin_merged, file = "cir_fin_merged.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(mrna_fin_merged , file = "mrna_fin_merged.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# keep higher than 400
cir_fin_merged_400set <- cir_fin_merged %>%
  filter(score > 400)

mrna_fin_merged_400set <- mrna_fin_merged %>%
  filter(score > 400)


# keep higher than 600
cir_fin_merged_600set <- cir_fin_merged %>%
  filter(score > 600)

mrna_fin_merged_600set <- mrna_fin_merged %>%
  filter(score > 600)

write.csv(cir_fin_merged_400set, file = "cir_fin_merged_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_fin_merged_600set, file = "cir_fin_merged_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(mrna_fin_merged_400set, file = "mrna_fin_merged_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(mrna_fin_merged_600set , file = "mrna_fin_merged_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

intersect(cir_fin_merged_400set$gene,mrna_fin_merged_400set$mRNA)
# 發現有3個gn 重複 在這邊+_cir > cir_fin_merged_400set2.csv
intersect(cir_fin_merged_600set$miRNA,mrna_fin_merged_600set$miRNA)
intersect(cir_fin_merged_400set$miRNA,mrna_fin_merged_400set$miRNA)

#保留cir id 有大於400 
#我要cir id 的read 根他interact mi 所以用merge家在一起
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2")
cir_fin=read.table(file = "cir_mi_cir_predict_fin.csv", header = TRUE, sep = ",")

cir_fin_400set <- cir_fin %>%
  filter(score > 400)
cir_fin_400set2= merge(x=cir_fin_400set, y=cir_att_list, by.x=c("gene"), by.y = c("host_gene"))

cir_fin_600set <- cir_fin %>%
  filter(score > 600)
cir_fin_600set2= merge(x=cir_fin_600set, y=cir_att_list, by.x=c("gene"), by.y = c("host_gene"))

#在手動分開作att list 根network就好
write.csv(cir_fin_400set2, file = "cir_fin_cir_id_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_fin_600set2 , file = "cir_fin_cir_id_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


#make mi att list>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mi_com= union(cir_fin_merged_400set$miRNA,mrna_fin_merged_400set$miRNA)
mi_com=as.data.frame(mi_com)
mi_com$mi_com= gsub( "-3p", "", mi_com$mi_com)
mi_com$mi_com= gsub( "-5p", "", mi_com$mi_com)
mi_com=unique(mi_com$mi_com)

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/pearson")

mi_pair2=read.table(file = "mi_pair2.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2")

mi_att_list=mi_pair2[ mi_pair2$mi_id %in% mi_com, ]
mi_com_not_in= mi_com[!mi_com %in% mi_pair2$mi_id]

#看完差多少enst 手動加回去 然後在去根fr拿 read count
write.csv(mi_att_list , file = "mi_att_list.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(mi_com_not_in , file = "mi_com_not_in.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#手動加入這些沒進去的mi gn 的enst 在來對照read count
#家完後 把一樣的enst 對應到把read 加進去
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/pearson")
mir_com_mrna_cir_fr2=read.table(file = "mir_com_mrna_cir_fr2.csv", header = TRUE, sep = ",")
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_list")
mi_att_list2=read.table(file = "mi_att_list.csv", header = TRUE, sep = ",")

mir_com_mrna_cir_fr2=mir_com_mrna_cir_fr2[,-c(2:9)]

mi_att_list3 <- merge(mi_att_list2, mir_com_mrna_cir_fr2, by.x = c("mi_enst"),by.y = c("enst"))

#家完後要把-3p -5p放回去 所以把他們對應到家回去
mi_com_400set= union(cir_fin_merged_400set$miRNA,mrna_fin_merged_400set$miRNA)
mi_com_400set=as.data.frame(mi_com400set)

write.csv(mi_com_400set , file = "mi_fin_com_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
#複製一個column 刪掉-3p -5p
mi_com_400set2=read.table(file = "mi_fin_com_400set.csv", header = TRUE, sep = ",")
mi_com_400set2$mi_remove= gsub( "-3p", "", mi_com_400set2$mi_remove)
mi_com_400set2$mi_remove= gsub( "-5p", "", mi_com_400set2$mi_remove)

mi_att_list4 <- merge(mi_com_400set2, mi_att_list3, by.x = c("mi_remove"), by.y = c("mi_id"), all.x = TRUE)

write.csv(mi_att_list4 , file = "mi_att_list_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

mi_com_600set=union(cir_fin_merged_600set$miRNA,mrna_fin_merged_600set$miRNA)
mi_att_list5= mi_att_list4[mi_att_list4$mi_com %in% mi_com_600set, ] 
write.csv(mi_att_list5 , file = "mi_att_list_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#make mi att list>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#make cir att list>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 我從這裏cir_com_mrna_mi_fr2 手動作出att_list_cir.csv
# 用 cir host gn 對應 400set根attlist 放上read count
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/group_score_subgn")
cir_fin_merged_400set=read.table(file = "cir_fin_merged_400set.csv", header = TRUE, sep = ",")
cir_fin_merged_600set=read.table(file = "cir_fin_merged_600set.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_list")
cir_att_list=read.table(file = "att_list_cir.csv", header = TRUE, sep = ",")

cir_com_400set=intersect(cir_fin_merged_400set$gene,cir_att_list$host_gene)
cir_att_list2= merge(x=cir_fin_merged_400set, y=cir_att_list, by.x=c("gene"), by.y = c("host_gene"))
cir_att_list2=cir_att_list2[,-c(2,3,5,6,10)]

#cir_att_list_cir_id= cir_att_list2[!duplicated(cir_att_list2$circRNA), ]
cir_att_list_cir_gn= cir_att_list2[!duplicated(cir_att_list2$gene), ]

#write.csv(cir_att_list_cir_id , file = "att_list_cir_id_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_att_list_cir_gn , file = "att_list_cir_gn_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

cir_com_600set=intersect(cir_fin_merged_600set$gene,cir_att_list$host_gene)
#cir_att_list_cir_id2= cir_att_list_cir_id[cir_att_list_cir_id$gene %in% cir_com_600set, ]  
cir_att_list_cir_gn2= cir_att_list_cir_gn[cir_att_list_cir_gn$gene %in% cir_com_600set, ]
# 之後在把host gn刪掉 cir_id用host gn 作對應
# 這邊cir id已經被merge過只剩一組 是錯的
#write.csv(cir_att_list_cir_id2 , file = "att_list_cir_id_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_att_list_cir_gn2 , file = "att_list_cir_gn_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


#make cir id att list  在上上面 以取得所有cir mi predict 400 600set
# 在重新作 cir id att list
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/group_score_subgn")
cir_fin_400set2 =read.table(file = "cir_fin_cir_id_400set.csv", header = TRUE, sep = ",")
cir_fin_600set2 =read.table(file = "cir_fin_cir_id_600set.csv", header = TRUE, sep = ",")

cir_att_list
cir_com_400set=intersect(cir_fin_400set2$gene,cir_att_list$host_gene)
cir_fin_400set3= merge(x=cir_fin_400set2, y=cir_att_list, by.x=c("gene"), by.y = c("host_gene"))

cir_com_600set=intersect(cir_fin_600set2$gene,cir_att_list$host_gene)
cir_fin_600set3= merge(x=cir_fin_600set2, y=cir_att_list, by.x=c("gene"), by.y = c("host_gene"))

write.csv(cir_fin_400set3 , file = "cir_fin_cir_id_400set2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_fin_600set3 , file = "cir_fin_cir_id_600set2.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#刪掉不要的col
#不要重複的cir id
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/group_score_subgn")
cir_fin_400set2 =read.table(file = "cir_fin_cir_id_400set2.csv", header = TRUE, sep = ",")
cir_fin_600set2 =read.table(file = "cir_fin_cir_id_600set2.csv", header = TRUE, sep = ",")

cir_att_list_cir_id= cir_fin_400set2[!duplicated(cir_fin_400set2$circRNA), ]
cir_att_list_cir_id2= cir_fin_600set2[!duplicated(cir_fin_600set2$circRNA), ]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_list")
write.csv(cir_att_list_cir_id , file = "att_list_cir_id_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_att_list_cir_id2 , file = "att_list_cir_id_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#make cir att list>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#make mrna att list>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#att_list_mrna.csv 是從 mrna_com_cir_mi_fr2 取出來

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/group_score_subgn")
mrna_fin_merged_400set=read.table(file = "mrna_fin_merged_400set.csv", header = TRUE, sep = ",")
mrna_fin_merged_600set=read.table(file = "mrna_fin_merged_600set.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_list")
mrna_att_list=read.table(file = "att_list_mrna.csv", header = TRUE, sep = "\t")

mrna_com_400set=intersect(mrna_fin_merged_400set$mRNA,mrna_att_list$gene_name)
mrna_com_600set=intersect(mrna_fin_merged_600set$mRNA,mrna_att_list$gene_name)

mrna_att_list2=mrna_att_list[mrna_att_list$gene_name %in% mrna_fin_merged_400set$mRNA, ]
mrna_att_list3=mrna_att_list[mrna_att_list$gene_name %in% mrna_fin_merged_600set$mRNA, ]

write.csv(mrna_att_list2 , file = "att_list_mrna_400set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(mrna_att_list3 , file = "att_list_mrna_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#make mrna att list>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#作 1010 2023 up dw gene regulate 在att list  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#在network做完 gsea 回來這邊作regulate table
# mrna cir  up dw gene regulate
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")

#把 network用到的att拿來用
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network")
att=read.table(file = "att_cir_gn_600set3_mrna_mi.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
upgn_mrna_gnname= read.table(file = "upgn_mrna_gnname.csv", header = TRUE, sep = ",")
dwgn_mrna_gnname= read.table(file = "dwgn_mrna_gnname.csv", header = TRUE, sep = ",")
upgn_cir_gnname= read.table(file = "upgn_cir_gnname.csv", header = TRUE, sep = ",")
dwgn_cir_gnname= read.table(file = "dwgn_cir_gnname.csv", header = TRUE, sep = ",")

att_up_mrna= att[att$gene %in% upgn_mrna_gnname$hgnc_symbol,]
att_dw_mrna= att[att$gene %in% dwgn_mrna_gnname$hgnc_symbol,]
att_up_cir= att[att$gene %in% upgn_cir_gnname$hgnc_symbol,]
att_dw_cir= att[att$gene %in% dwgn_cir_gnname$hgnc_symbol,]


att_up_mrna= mutate(att_up_mrna, regulated='up') 
att_dw_mrna= mutate(att_dw_mrna, regulated='dw') 
att_up_cir= mutate(att_up_cir, regulated='up') 
att_dw_cir= mutate(att_dw_cir, regulated='dw') 

write.csv(att_up_mrna , file = "att_up_mrna_regulate.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(att_dw_mrna , file = "att_dw_mrna_regulate.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(att_up_cir , file = "att_up_cir_regulate.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(att_dw_cir , file = "att_dw_cir_regulate.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#作 mi up dw gene regulate att list
setwd("/home/yuanpeng/yuan/R sarsCoV2/mirna_sequen/redo0920")

upgn_mi_gnname= read.table(file = "mirbase_mi_up.csv", header = TRUE, sep = ",")
dwgn_mi_gnname= read.table(file = "mirbase_mi_dw.csv", header = TRUE, sep = ",")

upgn_mi_gnname= grepl("^hsa" ,upgn_mi_gnname$mirbase_id)
dwgn_mi_gnname= grepl("^hsa" ,dwgn_mi_gnname$mirbase_id)


#作 up dw gene regulate 在att list  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>







#用手動把所有regulate table 加起來
#存在 setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
#存成att_updw_regulate_600set.csv


# 1010 2023 看共同mrna cir host gn >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 拿原本資料 total list 看 mrna cir共同出現的host gene是啥
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2")
mrna_mi_cir_total_list=read.table(file = "mrna_mi_cir_total_list.txt", header = FALSE, sep = ",")

colnames(mrna_mi_cir_total_list)= c("mRNA","miRNA","cir")
com_mrna_ci_total_list= data.frame()

for (i in 1:length(mrna_mi_cir_total_list$mRNA)) {
  if (mrna_mi_cir_total_list$mRNA[i] == mrna_mi_cir_total_list$cir[i]) {
    com_mrna_ci_total_list= rbind(com_mrna_ci_total_list,c(mrna_mi_cir_total_list$mRNA[i],mrna_mi_cir_total_list$cir[i]))
  }
}
colnames(com_mrna_ci_total_list)=c("mrna","cir")
freq=table(com_mrna_ci_total_list$mrna)
freq= freq[order(freq,decreasing = TRUE)]
barplot(freq)

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network/freq_com_host_gn")
write.csv(freq , file = "freq_m_cir_hostgn_total_list.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# 1010 2023 看共同mrna cir host gn >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# 1010 2023 看mrna-mi   cir-ni  highly conntectd node 100top >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

mrna_mi_node= paste0(mrna_mi_cir_total_list$mRNA,'_',mrna_mi_cir_total_list$miRNA) 
mrna_mi_node=as.data.frame(mrna_mi_node)
freq_mrna_mi_node= table(mrna_mi_node)
freq_mrna_mi_node= as.data.frame(freq_mrna_mi_node)

freq_mrna_mi_node2=freq_mrna_mi_node[order(freq_mrna_mi_node$Freq, decreasing = TRUE),]

cir_mi_node= paste0(mrna_mi_cir_total_list$miRNA,'_',mrna_mi_cir_total_list$cir)
cir_mi_node= as.data.frame(cir_mi_node)
freq_cir_mi_node= table(cir_mi_node)
freq_cir_mi_node= as.data.frame(freq_cir_mi_node)
freq_cir_mi_node2=freq_cir_mi_node[order(freq_cir_mi_node$Freq, decreasing = TRUE),]

library(tidyr)
freq_mrna_mi_node3= freq_mrna_mi_node2 %>%
  separate(
    col = mrna_mi_node,
    sep = "_",
    into = c("mRNA", "miRNA")
  )

freq_cir_mi_node3= freq_cir_mi_node2 %>%
  separate(
  col = cir_mi_node,
  sep = "_",
  into = c("miRNA", "circRNA")
  )

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network/freq_com_host_gn")
write.csv(freq_mrna_mi_node3 , file = "freq_mrna_mi_node.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(freq_cir_mi_node3 , file = "freq_cir_mi_node.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#各個rna 對所有connection  10/13
freq_mrna= table(mrna_mi_cir_total_list$mRNA)
freq_mrna= as.data.frame(freq_mrna)
freq_cir= table(mrna_mi_cir_total_list$cir)
freq_cir= as.data.frame(freq_cir)
freq_mi= table(mrna_mi_cir_total_list$miRNA)
freq_mi= as.data.frame(freq_mi)


freq_mrna2=freq_mrna[order(-freq_mrna$Freq),]
freq_mrna_100top= freq_mrna2[1:100,]
freq_cir2= freq_cir[order(-freq_cir$Freq),]
freq_cir_100top= freq_cir2[1:100,]
freq_mi2= freq_mi[order(-freq_mi$Freq), ]
freq_mi_100top=freq_mi2[1:100,]

write.csv(freq_mrna , file = "freq_mrna.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(freq_cir , file = "freq_cir.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(freq_mi , file = "freq_mi.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

write.csv(freq_mrna_100top , file = "freq_mrna_100top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(freq_cir_100top , file = "freq_cir_100top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(freq_mi_100top , file = "freq_mi_100top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

mrna_mi_= paste0(mrna_mi_cir_total_list$mRNA,'_',mrna_mi_cir_total_list$miRNA)


# 1010 2023 看mrna-mi   cir-mi  highly connected node >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



# 1013  取100top 共同score network list att list >>>>>>>>>>>>>>>

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network/freq_com_host_gn")
freq_mrna_100top=read.table(file = "freq_mrna_100top.csv", header = TRUE, sep = ",")
freq_cir_100top=read.table(file = "freq_cir_100top.csv", header = TRUE, sep = ",")
freq_mi_100top=read.table(file = "freq_mi_100top.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/group_score_subgn")
mrna_fin_merged=read.table(file = "mrna_fin_merged.csv", header = TRUE, sep = ",")
cir_fin_merged=read.table(file = "cir_fin_merged.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
mrna_fin_merged_100top= mrna_fin_merged[mrna_fin_merged$mRNA %in% freq_mrna_100top$Var1,]
cir_fin_merged_100top= cir_fin_merged[cir_fin_merged$gene %in% freq_cir_100top$Var1,]
#先filter mrna cir 在所有table中 在filter mi
mrna_fin_merged_100top= mrna_fin_merged_100top[mrna_fin_merged_100top$miRNA %in% freq_mi_100top$Var1,]
cir_fin_merged_100top= cir_fin_merged_100top[cir_fin_merged_100top$miRNA %in% freq_mi_100top$Var1,]

mrna_fin_merged_100top_600set= mrna_fin_merged_100top %>%
  filter(score > 600)

cir_fin_merged_100top_600set= cir_fin_merged_100top %>%
  filter(score > 600)

# 取得uni mi from mrna cir 的組合 後看共同mi 在作網路圖

mrna_cir_com_mi_100top_600set=intersect( mrna_fin_merged_100top_600set$miRNA,cir_fin_merged_100top_600set$miRNA)
mrna_cir_com_mi_100top_600set= unique(mrna_cir_com_mi_100top_600set)

mrna_fin_merged_100top_600set2=mrna_fin_merged_100top_600set[mrna_fin_merged_100top_600set$miRNA %in% mrna_cir_com_mi_100top_600set,]
cir_fin_merged_100top_600set2= cir_fin_merged_100top_600set[cir_fin_merged_100top_600set$miRNA %in% mrna_cir_com_mi_100top_600set,]

write.csv(mrna_fin_merged_100top_600set2 , file = "mrna_mi_100top_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_fin_merged_100top_600set2 , file = "cir_mi_100top_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
att_updw_regulate_600set= read.table(file = "att_updw_regulate_600set.csv", header = TRUE, sep = ",")

att_updw_regulate_600set_mrna=att_updw_regulate_600set[att_updw_regulate_600set$gene %in% mrna_fin_merged_100top_600set2$mRNA,]
att_updw_regulate_600set_cir=att_updw_regulate_600set[att_updw_regulate_600set$gene %in% cir_fin_merged_100top_600set2$gene,]
att_updw_regulate_600set_mi= att_updw_regulate_600set[att_updw_regulate_600set$gene %in% cir_fin_merged_100top_600set2$miRNA,]
#這邊確定有多少數字
length(intersect(att_updw_regulate_600set$gene,mrna_fin_merged_100top_600set2$mRNA))
length(intersect(att_updw_regulate_600set$gene,cir_fin_merged_100top_600set2$gene))
length(unique(mrna_fin_merged_100top_600set2$miRNA))
length(unique(cir_fin_merged_100top_600set2$miRNA))

att_updw_regulate_600set2= dplyr::bind_rows(att_updw_regulate_600set_cir,att_updw_regulate_600set_mrna, att_updw_regulate_600set_mi)

#存att list 手動把2組 mrna-mi  cir-mi table 家在一起
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
write.csv(att_updw_regulate_600set2 , file = "att_updw_regulate_100top_600set.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#作 mi


# 10 27 2023  >>>>>>>>>>>>>>>>>>>>
# 作top20 
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network/freq_com_host_gn")
freq_mrna_100top=read.table(file = "freq_mrna_100top.csv", header = TRUE, sep = ",")
freq_cir_100top=read.table(file = "freq_cir_100top.csv", header = TRUE, sep = ",")
freq_mi_100top=read.table(file = "freq_mi_100top.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top20")

freq_mrna_10top=freq_mrna_100top[1:10,]
freq_cir_10top= freq_cir_100top[1:10,]
freq_mi_10top= freq_mi_100top[1:10,]

mrna_fin_merged_10top= mrna_fin_merged[mrna_fin_merged$mRNA %in% freq_mrna_10top$Var1,]
cir_fin_merged_10top= cir_fin_merged[cir_fin_merged$gene %in% freq_cir_10top$Var1,]
#先filter mrna cir 在所有table中 在filter mi
mrna_fin_merged_10top= mrna_fin_merged_10top[mrna_fin_merged_10top$miRNA %in% freq_mi_10top$Var1,]
cir_fin_merged_10top= cir_fin_merged_10top[cir_fin_merged_10top$miRNA %in% freq_mi_10top$Var1,]

# 取得uni mi from mrna cir 的組合 後看共同mi 在作網路圖

mrna_cir_com_mi_10top=intersect( mrna_fin_merged_10top$miRNA,cir_fin_merged_10top$miRNA)
mrna_cir_com_mi_10top= unique(mrna_cir_com_mi_10top)

mrna_fin_merged_10top2=mrna_fin_merged_10top[mrna_fin_merged_10top$miRNA %in% mrna_cir_com_mi_10top,]
cir_fin_merged_10top2= cir_fin_merged_10top[cir_fin_merged_10top$miRNA %in% mrna_cir_com_mi_10top,]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top10")
write.csv(mrna_fin_merged_10top2 , file = "mrna_mi_10top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_fin_merged_10top2 , file = "cir_mi_10topset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# att list 拿就的來用 發現有少一些 手動家進去
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
att_updw_regulate_600set= read.table(file = "att_updw_regulate_600set.csv", header = TRUE, sep = ",")
length(intersect(att_updw_regulate_600set$gene,mrna_fin_merged_20top2$mRNA)) 
length(intersect(att_updw_regulate_600set$gene,cir_fin_merged_20top2$gene)) 
length(intersect(att_updw_regulate_600set$gene,cir_fin_merged_20top2$miRNA))

att_updw_regulate_20set_mrna=att_updw_regulate_600set[att_updw_regulate_600set$gene %in% mrna_fin_merged_20top2$mRNA,]
att_updw_regulate_20set_cir=att_updw_regulate_600set[att_updw_regulate_600set$gene %in% cir_fin_merged_20top2$gene,]
att_updw_regulate_20set_mi= att_updw_regulate_600set[att_updw_regulate_600set$gene %in% cir_fin_merged_20top2$miRNA,]

att_updw_regulate_20set2= dplyr::bind_rows(att_updw_regulate_20set_cir,att_updw_regulate_20set_mrna, att_updw_regulate_20set_mi)

#存att list 手動把2組 mrna-mi  cir-mi table 家在一起
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top10")
write.csv(att_updw_regulate_20set2 , file = "att_updw_regulate_10top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


#top 20 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

freq_mrna_20top=freq_mrna_100top[1:20,]
freq_cir_20top= freq_cir_100top[1:20,]
freq_mi_20top= freq_mi_100top[1:20,]

mrna_fin_merged_20top= mrna_fin_merged[mrna_fin_merged$mRNA %in% freq_mrna_20top$Var1,]
cir_fin_merged_20top= cir_fin_merged[cir_fin_merged$gene %in% freq_cir_20top$Var1,]
#先filter mrna cir 在所有table中 在filter mi
mrna_fin_merged_20top= mrna_fin_merged_20top[mrna_fin_merged_20top$miRNA %in% freq_mi_20top$Var1,]
cir_fin_merged_20top= cir_fin_merged_20top[cir_fin_merged_20top$miRNA %in% freq_mi_20top$Var1,]

# 取得uni mi from mrna cir 的組合 後看共同mi 在作網路圖

mrna_cir_com_mi_20top=intersect( mrna_fin_merged_20top$miRNA,cir_fin_merged_20top$miRNA)
mrna_cir_com_mi_20top= unique(mrna_cir_com_mi_20top)

mrna_fin_merged_20top2=mrna_fin_merged_20top[mrna_fin_merged_20top$miRNA %in% mrna_cir_com_mi_20top,]
cir_fin_merged_20top2= cir_fin_merged_20top[cir_fin_merged_20top$miRNA %in% mrna_cir_com_mi_20top,]

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top20")
write.csv(mrna_fin_merged_20top2 , file = "mrna_mi_20top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(cir_fin_merged_20top2 , file = "cir_mi_20top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# att list 拿就的來用 發現有少一些 手動家進去 

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
att_updw_regulate_600set= read.table(file = "att_updw_regulate_600set.csv", header = TRUE, sep = ",")
length(intersect(att_updw_regulate_600set$gene,mrna_fin_merged_20top2$mRNA)) 
length(intersect(att_updw_regulate_600set$gene,cir_fin_merged_20top2$gene)) 
length(intersect(att_updw_regulate_600set$gene,cir_fin_merged_20top2$miRNA))

att_updw_regulate_20set_mrna=att_updw_regulate_600set[att_updw_regulate_600set$gene %in% mrna_fin_merged_20top2$mRNA,]
att_updw_regulate_20set_cir=att_updw_regulate_600set[att_updw_regulate_600set$gene %in% cir_fin_merged_20top2$gene,]
att_updw_regulate_20set_mi= att_updw_regulate_600set[att_updw_regulate_600set$gene %in% cir_fin_merged_20top2$miRNA,]

att_updw_regulate_20set2= dplyr::bind_rows(att_updw_regulate_20set_cir,att_updw_regulate_20set_mrna, att_updw_regulate_20set_mi)

#存att list 手動把2組 mrna-mi  cir-mi table 家在一起
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top20")
write.csv(att_updw_regulate_20set2 , file = "att_updw_regulate_20top.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)





