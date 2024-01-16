library(ggplot2)
library(ggraph)
library(tidygraph)
library(igraph)
library(dbplyr)
library(tidyr)
library(caret)

setwd("/home/yuanpeng/yuan/mirnada/newpredict")

df=read.table(file = "test2.csv", header = TRUE, sep = "\t")
att= read.table(file = "att_list.csv", header = TRUE, sep = ",")


ci=ci[,c('Seq1','logtotalread')] 
mi=mi[,c('Seq1','logtotalread')] 

#作table 選擇 target source socre
df.el=df%>%select(Seq1,Seq2,Tot.Score)%>%
  group_by(Seq1,Seq2)%>%
  expand(edge=c(1:Tot.Score))%>%select(-edge)
head(df.el)

#作matrix 給畫圖
got1=graph_from_edgelist(df.el%>%as.matrix(),directed = FALSE)
#簡單畫圖 要跑漫九的
plot(got1,edge.arrow.size=.5, vertex.color="gold", vertex.size=3, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.2)
#>>>>>>>>>>
#這邊算totscore 把他全部變成1 然後看出現幾次1 算sum畫圖 但我要保留edge的寬度
E(got1)$Tot.Score=rep(1,ecount(got1))
got1s <- igraph::simplify( got1, remove.multiple = T, remove.loops = F, 
                           edge.attr.comb=c(Tot.Score="sum"))

#>>>>>>>>>>>>>>>
got1s <- igraph::simplify(got1, remove.multiple = TRUE, remove.loops = FALSE)
plot(got1s,edge.arrow.size=.5, vertex.color="gold", vertex.size=3, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=.5, vertex.label.dist=2, edge.curved=0.5,layout=layout_with_lgl)

got2 = make_graph( edges = df[, c("Seq1", "Seq2")],
                   directed = FALSE)


got2 <- as_tbl_graph(got1s)

got2 %>% create_notable('zachary') %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree())

create_notable('zachary') %>%
  activate(nodes) %>% 
  mutate(group = as.factor(group_infomap())) %>% # Creates a `group` variable based on the infomap algorithm
  ggraph(layout = 'stress') +
  geom_edge_fan(width = .2, color = 'lightblue') + 
  geom_node_point(aes(color = group)) + 
  coord_fixed() + 
  theme_graph()


relation= data.frame( from= df$Seq1, to= df$Seq2, edge_score=df$Tot.Score )

got3 = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

got3 %>%
  ggraph(layout = 'stress') +
  geom_edge_fan(width = .5, color = 'gray') +
  geom_node_point(aes(color = , size = 3)) +
  scale_color_viridis_c()

node_data(got3)
got3


graph <- graph_from_data_frame(tbl_links,tbl_nodes, directed = TRUE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$logtotalread))
got3 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = edge_score)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount ), size = 0.5) +
  geom_node_text(aes(label = att$Seq1, colour = group), repel = TRUE, size = 0.5) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

# network_m_cir_mi2.csv 把所有重複都丟進去跑 這裏把部重複的丟進去跑 反正多的都不會在圖裏

df=read.table(file = "network_m_cir_mi.csv", header = TRUE, sep = ",")
att= read.table(file = "att_list.csv", header = TRUE, sep = ",")

relation2= data.frame( from= df$Seq1, to= df$Seq2, edge_score=df$logtscore )

got4 = graph_from_data_frame(relation2, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount2=as.integer(as.factor(att$logtotalread))
node_degree <- degree(got4, mode="all")

got4 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = edge_score)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount2 ), size = 5) +
  geom_node_text(aes(label = att$Seq1, colour = group), repel = TRUE, size = 5) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

logreadcount3=as.integer(att$logtotalread)
got4 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = edge_score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 , size = node_degree)) +
  geom_node_text(aes(label = att$Seq1, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色


# >>>>>>>> 重跑 所有正常id  把沒有的cir刪掉 >>>>>>>>>>>>>>>>>>>

df=read.table(file = "network_m_cir_mi3.csv", header = TRUE, sep = ",")
att= read.table(file = "att_list2.csv", header = TRUE, sep = ",")

relation2= data.frame( from= df$Seq1, to= df$Seq2, edge_score=df$logscore )

got4 = graph_from_data_frame(relation2, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount2=as.integer(as.factor(att$logtotalread))
node_degree <- degree(got4, mode="all")

got4 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = edge_score)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount2 ), size = 5) +
  geom_node_text(aes(label = att$Seq1, colour = group), repel = TRUE, size = 5) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

logreadcount3=as.integer(att$logtotalread)
got4 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = edge_score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 , size = node_degree)) +
  geom_node_text(aes(label = att$Seq1, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色


# >>>>>>>>>>>> group 把所有cir加進同一組host gn >>>>>>>>>>>>>>

cir_net=read.table(file = "test3_cir.csv", header = TRUE, sep = "\t")
mrna_net=read.table(file = "test2_mrna.csv", header = TRUE, sep = "\t")

cir_net2 <- cir_net %>% group_by(Seq1,Seq2_gnname2 ) %>% 
  summarise(Tot.Score=sum(Tot.Score))


test_cir <- cir_net %>%
  group_by(Seq1, Seq2_gnname2) %>%
  count() %>%
  rename(node_degree = n)

mrna_net2 <- mrna_net %>% group_by(Seq1,Seq2 ) %>% 
  summarise(Tot.Score=sum(Tot.Score))

test_mrna <- mrna_net %>%
  group_by(Seq1, Seq2) %>%
  count() %>%
  rename(node_degree = n)

log_cir=as.data.frame(log(cir_net2$Tot.Score))
log_mrna=as.data.frame(log(mrna_net2$Tot.Score))

cir_net3= cbind(cir_net2,log_cir)
mrna_net3= cbind(mrna_net2,log_mrna)

test_cir$Seq2_gnname2= paste0(test_cir$Seq2_gnname2, '_circRNA')
test_mrna$Seq2= paste0(test_mrna$Seq2, '_mRNA')


write.csv(test_cir, file = "cir_combine_node.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(test_mrna, file = "mrna_combine_node.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#>>>>>>>>>>>>>>>>>>>>>make att list node >>>>>>>>>>>>>>>>

att= read.table(file = "att_list.csv", header = TRUE, sep = ",")
cir_net=read.table(file = "test3_cir.csv", header = TRUE, sep = "\t")

new_att=data.frame()
#取一樣cir id 把cir id放geneid
for (i in 1:length(att$Seq1)) {
  for (q in 1:length(cir_net$Seq2)) {
    if (att$Seq1[i] == cir_net$Seq2[q]) {
      new_att= rbind(new_att,c(att$Seq1[i],cir_net$Seq2_gnname2[q] ,att$group[i],att$logtotalread[i])) 
    }
  }
}

colnames(new_att)= c("cir_id","host_gene",'group','logtotalread')

write.csv(new_att, file = "att_list_cir2.csv", sep=",",eol = "\n", row.names=T, col.names= TRUE)

#>>>>>>>>>>>>>>>> graph combine node >>>>>>>>>>>>>>

# network imporved
df2=read.table(file = "network_m_cir_mi_combine_node.csv", header = TRUE, sep = ",")
att2=read.table(file = "att_list_combine_node.csv", header = TRUE, sep = ",")

relation3= data.frame( from= df2$Seq1, to= df2$Seq2_gnname, sub_gene=df2$sub_gene )

got5 = graph_from_data_frame(relation3, att2,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount2=as.integer(as.factor(att2$logtotalread))
node_degree <- degree(got5, mode="all")


got5 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 ), size = 3) +
  geom_node_text(aes(label = att2$host_gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

logreadcount3=as.integer(att2$logtotalread)
got5 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 , size = node_degree)) +
  geom_node_text(aes(label = att2$host_gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mrna_mi_predict=read.table(file = "mrna_mi_predict55.csv", header = TRUE, sep = "\t")
cir_mi_predict=read.table(file = "test3_cir_mi_predict.csv", header = TRUE, sep = "\t")

mrna_mi_predict2=data.frame()
cir_mi_predict2=data_frame()

com_gn_name=intersect(mrna_mi_predict$Seq2,cir_mi_predict$Seq2_gnname2)
com_gn_name= as.data.frame(com_gn_name)
cir_mi_predict2=cir_mi_predict[cir_mi_predict$Seq2_gnname2 %in% com_gn_name$`intersect(mrna_mi_predict$Seq2, cir_mi_predict$Seq2_gnname2)`,]
mrna_mi_predict2=mrna_mi_predict[mrna_mi_predict$Seq2 %in% com_gn_name$`intersect(mrna_mi_predict$Seq2, cir_mi_predict$Seq2_gnname2)`,]



test_cir <- cir_mi_predict2 %>%
  group_by(Seq1, Seq2_gnname2) %>%
  count() %>%
  rename(node_degree = n)


test_mrna <- mrna_mi_predict2 %>%
  group_by(Seq1, Seq2) %>%
  count() %>%
  rename(node_degree = n)

test_cir$Seq2_gnname2= paste0(test_cir$Seq2_gnname2, '_circRNA')
test_mrna$Seq2= paste0(test_mrna$Seq2, '_mRNA')


write.csv(test_cir, file = "cir_combine_node_comm_mrna.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(test_mrna, file = "mrna_combine_node_comm_cir.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#做完共同network 作att list
att2=read.table(file = "att_list_combine_node.csv", header = TRUE, sep = "\t")

test_mrna2= unique(test_mrna$Seq2)
test_cir2= unique(test_cir$Seq2_gnname2)
com_gn_name2= c(test_mrna2, test_cir2)

att3= att2[att2$host_gene %in% com_gn_name2,]

write.csv(att3, file = "att_list_combine_node_common.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

#做好 network att list 刪掉沒有的 作network圖
df3=read.table(file = "network_m_cir_mi_combine_node_common.csv", header = TRUE, sep = ",")
att3=read.table(file = "att_list_combine_node_common.csv", header = TRUE, sep = ",")

relation4= data.frame( from= df3$Seq1, to= df3$Seq2, sub_gene=df3$sub_gene )

got6 = graph_from_data_frame(relation4, att3,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount3=as.integer(as.factor(att3$logtotalread))
node_degree <- degree(got6, mode="all")


got6 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 ), size = 3) +
  geom_node_text(aes(label = att3$host_gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

logreadcount3=as.integer(att3$logtotalread)
got6 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 , size = node_degree)) +
  geom_node_text(aes(label = att3$host_gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

# delete unique node
df3=read.table(file = "network_m_cir_mi_combine_node_common2.csv", header = TRUE, sep = ",")
att3=read.table(file = "att_list_combine_node_common2.csv", header = TRUE, sep = ",")

relation4= data.frame( from= df3$Seq1, to= df3$Seq2, sub_gene=df3$sub_gene )

got6 = graph_from_data_frame(relation4, att3,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount3=as.integer(as.factor(att3$logtotalread))
node_degree <- degree(got6, mode="all")


got6 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 ), size = 3) +
  geom_node_text(aes(label = att3$host_gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

logreadcount3=as.integer(att3$logtotalread)
got6 %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount3 , size = node_degree)) +
  geom_node_text(aes(label = att3$host_gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色


#>>>>>>>>>>>>>>>>>>>10/05  new network >>>>>>>>>>>>>>>>>>>>>>
#cir combine node  cir gn group
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network")

df=read.table(file = "network_combine_node_cir_gn_400set.csv", header = TRUE, sep = ",")
att=read.table(file = "att_list_cir_gn_400set.csv", header = TRUE, sep = ",")
#sub_gene =基因數量
relation= data.frame( from= df$seq1, to= df$seq2, sub_gene=df$node_degree )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")


got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount ), size = 3) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

# too messy

# try 600set

df=read.table(file = "network_combine_node_cir_gn_600set.csv", header = TRUE, sep = "\t")
att=read.table(file = "att_list_cir_gn_600set.csv", header = TRUE, sep = ",")
#sub_gene =基因數量
relation= data.frame( from= df$seq1, to= df$seq2, sub_gene=df$node_degree )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount ), size = 3) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色
#save network_cir_gn_600set
#remove unique pair  把多出來獨立的pair山了

df=read.table(file = "network_combine_node_cir_gn_600set2.csv", header = TRUE, sep = ",")
att=read.table(file = "att_list_cir_gn_600set2.csv", header = TRUE, sep = ",")
#sub_gene =基因數量
relation= data.frame( from= df$seq1, to= df$seq2, sub_gene=df$node_degree )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount ), size = 3) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色
# network_cir_gn_600set2
# remove mirna only showed once 在把只有出現1次的mi刪掉

df=read.table(file = "network_combine_node_cir_gn_600set3.csv", header = TRUE, sep = ",")
att=read.table(file = "att_list_cir_gn_600set3.csv", header = TRUE, sep = ",")
#sub_gene =基因數量
relation= data.frame( from= df$seq1, to= df$seq2, sub_gene=df$node_degree )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount ), size = 3) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = sub_gene),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色
# save network_cir_gn_600set3



#>>>>>>>>>>>>>>>>>>>10/05  new network >>>>>>>>>>>>>>>>>>>>>>
#cir combine node  cir gn group




#>>>>>>>>>>>>>>>>>>>10/05  new network >>>>>>>>>>>>>>>>>>>>>>
#cir   cir id group 600set

df=read.table(file = "network_cir_id_600set.csv", header = TRUE, sep = "\t")
att=read.table(file = "att_list_cir_id_600set.csv", header = TRUE, sep = ",")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount ), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色
#save network_cir_id_600set.pdf
#刪掉單獨的pair

df=read.table(file = "network_cir_id_600set2.csv", header = TRUE, sep = ",")
att=read.table(file = "att_list_cir_id_600set2.csv", header = TRUE, sep = ",")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group , alpha=logreadcount), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白

#刪掉出現一次的mi
#我手動在excel 看只出現1次的mi 刪掉 作一次圖 在把獨立的mi從att刪掉
df=read.table(file = "network_cir_id_600set3.csv", header = TRUE, sep = ",")
#att=read.table(file = "att_list_cir_id_600set2.csv", header = TRUE, sep = ",")
att=read.table(file = "att_list_cir_id_600set3.csv", header = TRUE, sep = ",")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("red", "blue", "green")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group , alpha=logreadcount), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色
# save network_cir_id_600set.pdf

#>>>>>>>>>>>>>>>>>>>10/05  new network >>>>>>>>>>>>>>>>>>>>>>
#cir   cir id group


#只保留有根cir 互動的mi 在network李 所以mrna數量變少
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network")

df=read.table(file = "network_cir_gn_600set32.csv", header = TRUE, sep = ",")
att=read.table(file = "att_cir_gn_600set3_mrna_mi.csv", header = TRUE, sep = ",")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("#bb5548", "#745399", "#715c1f")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group , alpha=logreadcount), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

# save network_cir_gn_600set32.pdf



# 10/10 2023  作up dw regulate 的 network圖
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network")
df=read.table(file = "network_cir_gn_600set32.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")
att=read.table(file = "att_updw_regulate_600set.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network/regulate")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("#c97586", "#86C975", "#7586C9")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group , alpha=logreadcount,shape = regulated), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色


# 把node degree  normalize 用preprocess 要下載package, 設定範圍1～20
# 加上 size = norm_scale$node_degree, shape = regulated
# shape 有圖是可以選 scale_shape_manual(values = c(15, 17, 19)) 這編選 第15 17 19的圖案
# 因為顏色太透明 所以用 scale_alpha_continuous(range = c(0.3, 1)) 設定顏色 
# scale_alpha_manual 部之有沒有這個
library(caret)

group_colors <- c("#57A27A", "#A27A57", "#7A57A2")
group_colors <- c("#BE5A6E", "#6EBE5A", "#5A6EBE")


process <- preProcess(as.data.frame(node_degree), method=c("range"),rangeBounds = c(1, 20))
norm_scale <- predict(process, as.data.frame(node_degree))

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount, size = norm_scale$node_degree, shape = regulated)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  scale_alpha_continuous(range = c(0.3, 1))+
  labs(size = "node_degree", edge_alpha="miranda_score")+
  theme_void() #把背景便白色



# 10/10 2023 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# finish network then put mrna cirrna into GSEA
# ori 900 mrna gn
setwd("/home/yuanpeng/yuan/R sarsCoV2/result/redo_0920")
upgn_mrna= read.table(file = "upgn_fr_mrna.csv", header = TRUE, sep = ",")
dwgn_mrna=read.table(file = "dwgn_fr_mrna.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/R sarsCoV2/cirRNA_ensg/new2/redo_0920")
upgn_cir= read.table(file = "upgn_fr_cir.csv", header = TRUE, sep = ",")
dwgn_cir=read.table(file = "dwgn_fr_cir.csv", header = TRUE, sep = ",")

setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/GSEA")

upgn_mrna=upgn_mrna[,1]
dwgn_mrna=dwgn_mrna[,1]
upgn_cir=upgn_cir[,1]
dwgn_cir=dwgn_cir[,1]

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upgn_mrna_ID=biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name'),
                       filters = 'ensembl_gene_id', 
                       values = upgn_mrna,
                       mart = mart)

dwgn_mrna_ID=biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name'),
                            filters = 'ensembl_gene_id', 
                            values = dwgn_mrna,
                            mart = mart)

upgn_cir_ID=biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol','ensembl_transcript_id'),
                            filters = 'ensembl_transcript_id', 
                            values = upgn_cir,
                            mart = mart)

dwgn_cir_ID=biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol','ensembl_transcript_id'),
                            filters = 'ensembl_transcript_id', 
                            values = dwgn_cir,
                            mart = mart)

write.csv(upgn_mrna_ID, file = "upgn_mrna_gnname.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_mrna_ID, file = "dwgn_mrna_gnname.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(upgn_cir_ID, file = "upgn_cir_gnname.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)
write.csv(dwgn_cir_ID, file = "dwgn_cir_gnname.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

# gprofiler for GSEA

upgn_mrna_gnset= gprofiler2::gost(upgn_mrna_ID$hgnc_symbol, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
dwgn_mrna_gnset= gprofiler2::gost(dwgn_mrna_ID$hgnc_symbol, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
upgn_cir_gnset= gprofiler2::gost(upgn_cir_ID$hgnc_symbol, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)
dwgn_cir_gnset= gprofiler2::gost(dwgn_cir_ID$hgnc_symbol, correction_method="fdr", sources =c("GO:BP","GO:MF"),evcodes = TRUE)

upgn_mrna_gnset=upgn_mrna_gnset$result
dwgn_mrna_gnset=dwgn_mrna_gnset$result
upgn_cir_gnset=upgn_cir_gnset$result
dwgn_cir_gnset=dwgn_cir_gnset$result

# select only 10-1000 of term size
upgn_mrna_gnset2<- data.frame()
for (i in 1:length(upgn_mrna_gnset$term_size)){
  if(upgn_mrna_gnset$term_size[i]>10 && upgn_mrna_gnset$term_size[i]<1000){
    print(i)
    print(upgn_mrna_gnset$term_size[i])
    upgn_mrna_gnset2=rbind(upgn_mrna_gnset2,c(upgn_mrna_gnset$term_size[i],upgn_mrna_gnset$source[i],
                                              upgn_mrna_gnset$p_value[i],upgn_mrna_gnset$term_id[i]
                                              ,upgn_mrna_gnset$term_name[i],upgn_mrna_gnset$query_size[i]
                                              ,upgn_mrna_gnset$intersection_size[i],upgn_mrna_gnset$intersection[i]))
  }
}

colnames(upgn_mrna_gnset2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
# order by p value  order(as.numeric) let 0 e-4 become numeric and can be ordered
upgn_mrna_gnset2=upgn_mrna_gnset2[order(as.numeric(upgn_mrna_gnset2$p_value)),]
# change num size of p value to become only 5 number > 0.00001
#format 可以用科學記號
#upgn_mrna_gnset2$p_value=format(as.numeric(upgn_mrna_gnset2$p_value) , digits = 5, scientific = TRUE)
write.csv(upgn_mrna_gnset2, file = "upgn_mrna_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

dwgn_mrna_gnset2<- data.frame()
for (i in 1:length(dwgn_mrna_gnset$term_size)){
  if(dwgn_mrna_gnset$term_size[i]<1000 && dwgn_mrna_gnset$term_size[i]>10){
    print(i)
    print(dwgn_mrna_gnset$term_size[i])
    dwgn_mrna_gnset2=rbind(dwgn_mrna_gnset2,c(dwgn_mrna_gnset$term_size[i],dwgn_mrna_gnset$source[i],
                                              dwgn_mrna_gnset$p_value[i],dwgn_mrna_gnset$term_id[i]
                                              ,dwgn_mrna_gnset$term_name[i],dwgn_mrna_gnset$query_size[i]
                                              ,dwgn_mrna_gnset$intersection_size[i],dwgn_mrna_gnset$intersection[i]))
  }
}

colnames(dwgn_mrna_gnset2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
dwgn_mrna_gnset2=dwgn_mrna_gnset2[order(as.numeric(dwgn_mrna_gnset2$p_value)),]
# change num size of p value to become only 5 number > 0.00001
#dwgn_mrna_gnset2$p_value=round(as.numeric(dwgn_mrna_gnset2$p_value) , digits = 5)
write.csv(dwgn_mrna_gnset2, file = "dwgn_mrna_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

upgn_cir_gnset2<- data.frame()
for (i in 1:length(upgn_cir_gnset$term_size)){
  if(upgn_cir_gnset$term_size[i]>10 && upgn_cir_gnset$term_size[i]<1000){
    print(i)
    print(upgn_cir_gnset$term_size[i])
    upgn_cir_gnset2=rbind(upgn_cir_gnset2,c(upgn_cir_gnset$term_size[i],upgn_cir_gnset$source[i],
                                              upgn_cir_gnset$p_value[i],upgn_cir_gnset$term_id[i]
                                              ,upgn_cir_gnset$term_name[i],upgn_cir_gnset$query_size[i]
                                              ,upgn_cir_gnset$intersection_size[i],upgn_cir_gnset$intersection[i]))
  }
}

colnames(upgn_cir_gnset2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
# order by p value  order(as.numeric) let 0 e-4 become numeric and can be ordered
upgn_cir_gnset2=upgn_cir_gnset2[order(as.numeric(upgn_cir_gnset2$p_value)),]
# change num size of p value to become only 5 number > 0.00001
#upgn_cir_gnset2$p_value=format(as.numeric(upgn_cir_gnset2$p_value) , digits = 5, scientific = TRUE)
write.csv(upgn_cir_gnset2, file = "upgn_cir_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)

dwgn_cir_gnset2<- data.frame()
for (i in 1:length(dwgn_cir_gnset$term_size)){
  if(dwgn_cir_gnset$term_size[i]>10 && dwgn_cir_gnset$term_size[i]<1000){
    print(i)
    print(dwgn_cir_gnset$term_size[i])
    dwgn_cir_gnset2=rbind(dwgn_cir_gnset2,c(dwgn_cir_gnset$term_size[i],dwgn_cir_gnset$source[i],
                                            dwgn_cir_gnset$p_value[i],dwgn_cir_gnset$term_id[i]
                                            ,dwgn_cir_gnset$term_name[i],dwgn_cir_gnset$query_size[i]
                                            ,dwgn_cir_gnset$intersection_size[i],dwgn_cir_gnset$intersection[i]))
  }
}

colnames(dwgn_cir_gnset2 )=c('term_size','source',"p_value",'term_ID',"term_name","query_size","intersection_size","intersection")
# order by p value  order(as.numeric) let 0 e-4 become numeric and can be ordered
dwgn_cir_gnset2=dwgn_cir_gnset2[order(as.numeric(dwgn_cir_gnset2$p_value)),]
# change num size of p value to become only 5 number > 0.00001
#dwgn_cir_gnset2$p_value=format(as.numeric(dwgn_cir_gnset2$p_value) , digits = 5, scientific = TRUE)
write.csv(dwgn_cir_gnset2, file = "dwgn_cir_gnset.csv", sep=",",eol = "\n", row.names=F, col.names= TRUE)


# 10/17 2023  作up dw regulate 的 network圖  100top >>>>>>>>>>>>>>>>>>>>>>>
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top")
df=read.table(file = "network_gn_100top_600set.csv", header = TRUE, sep = ",")

att=read.table(file = "att_updw_regulate_100top_600set.csv", header = TRUE, sep = ",")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("#c97586", "#86C975", "#7586C9")  # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group , alpha=logreadcount,shape = regulated), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色


# 把node degree  normalize 用preprocess 要下載package, 設定範圍1～20
# 加上 size = norm_scale$node_degree, shape = regulated
# shape 有圖是可以選 scale_shape_manual(values = c(15, 17, 19)) 這編選 第15 17 19的圖案
# 因為顏色太透明 所以用 scale_alpha_continuous(range = c(0.3, 1)) 設定顏色 
# scale_alpha_manual 部之有沒有這個
library(caret)

group_colors <- c("#57A27A", "#A27A57", "#7A57A2")
group_colors <- c("#BE5A6E", "#6EBE5A", "#5A6EBE")
group_colors <- c("#4E1CC8", "#D63E29", "#96C81C")

process <- preProcess(as.data.frame(node_degree), method=c("range"),rangeBounds = c(1, 20))
norm_scale <- predict(process, as.data.frame(node_degree))

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount, size = norm_scale$node_degree, shape = regulated)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  scale_alpha_continuous(range = c(0.3, 1))+
  labs(size = "node_degree", edge_alpha="miranda_score")+
  theme_void() #把背景便白色 



# 10/27 2023  作up dw regulate 的 network圖  10top >>>>>>>>>>>>>>>>>>>>>>>
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top10")
df=read.table(file = "network_gn_10top.csv", header = TRUE, sep = ",")

att=read.table(file = "att_updw_regulate_10top.csv", header = TRUE, sep = ",")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()
group_colors <- c("#c97586", "#47885e", "#4c6cb3") 
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group , alpha=logreadcount,shape = regulated), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

process <- preProcess(as.data.frame(node_degree), method=c("range"),rangeBounds = c(1, 20))
norm_scale <- predict(process, as.data.frame(node_degree))
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount, size = norm_scale$node_degree, shape = regulated)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  scale_alpha_continuous(range = c(0.3, 1))+
  labs(size = "node_degree", edge_alpha="miranda_score")+
  theme_void() #把背景便白色 

# 1102 老師說用 top10 >>>>>>>>>>>>>>>>>>>>>>>
# 作有hsa 根梅 hsa的組
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount, size = node_degree, shape = regulated)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 5) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  scale_alpha_continuous(range = c(0.3, 1))+
  labs(size = "node_degree", edge_alpha="miranda_score")+
  theme_void(base_size = 13) #把背景便白色 

# 移除 hsa from mir name
df=read.table(file = "network_gn_10top_nohsa.csv", header = TRUE, sep = ",")
att=read.table(file = "att_updw_regulate_10top_nohsa.csv", header = TRUE, sep = ",")
relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )
got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount, size = node_degree, shape = regulated)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 5) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  scale_alpha_continuous(range = c(0.3, 1))+
  labs(size = "node_degree", edge_alpha="miranda_score")+
  theme_void(base_size = 13) #把背景便白色 

freq_mi= table(df$seq1)
freq_mi=freq_mi[order(freq_mi, decreasing = TRUE)]
barplot(freq_mi)

att_m= att[att$group == "mRNA", ]
att_c= att[att$group == "circRNA", ]
df_m= merge(x=df, y=att_m, by.x=c("seq2"), by.y = c("gene"))
df_c= merge(x=df, y=att_c, by.x=c("seq2"), by.y = c("gene"))

freq_m= table(df_m$seq2)
freq_m=freq_m[order(freq_m, decreasing = TRUE)]
barplot(freq_m)

freq_c= table(df_c$seq2)
freq_c=freq_c[order(freq_c, decreasing = TRUE)]
barplot(freq_c)


# 1102 老師說用 top10 >>>>>>>>>>>>>>>>>>>>>>>

# 10/27 2023  作up dw regulate 的 network圖  20top >>>>>>>>>>>>>>>>>>>>>>>
setwd("/home/yuanpeng/yuan/Rsarscov_a549_network2/network_100top/top20")
df=read.table(file = "network_gn_20top.csv", header = TRUE, sep = ",")

att=read.table(file = "att_updw_regulate_20top.csv", header = TRUE, sep = ",")

relation= data.frame( from= df$seq1, to= df$seq2, score=df$score )

got = graph_from_data_frame(relation, att,directed = FALSE) %>% as_tbl_graph()

group_colors <- c("#c97586", "#47885e", "#4c6cb3")   # Define your group colors here
logreadcount=as.integer(as.factor(att$log_total_read))
node_degree <- degree(got, mode="all")

got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score)) +
  geom_node_point(aes(colour = group, fill = group , alpha=logreadcount,shape = regulated), size = 3) +
  geom_node_text(aes(label = att$gene , colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白

logreadcount=as.integer(att$log_total_read)
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount , size = node_degree)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_colour_manual(values = group_colors) +  # Apply group colors
  theme_void() #把背景便白色

process <- preProcess(as.data.frame(node_degree), method=c("range"),rangeBounds = c(1, 20))
norm_scale <- predict(process, as.data.frame(node_degree))
got %>%
  ggraph(layout = 'stress') +
  geom_edge_link(aes(edge_alpha = score),force_flip = TRUE) +
  geom_node_point(aes(colour = group, fill = group, alpha=logreadcount, size = norm_scale$node_degree, shape = regulated)) +
  geom_node_text(aes(label = att$gene, colour = group), repel = TRUE, size = 3) +
  scale_shape_manual(values = c(15, 17, 19))+
  scale_colour_manual(values = group_colors) +  # Apply group colors
  scale_alpha_continuous(range = c(0.3, 1))+
  labs(size = "node_degree", edge_alpha="miranda_score")+
  theme_void() #把背景便白色 


