#1.数据读取----
#加载包
library(Hmisc)
library(dplyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(ComplexHeatmap)
#1.1读取注释表----
ko1_4 = read.table("E:/PCS/KO1-4.txt",header = TRUE, sep = "\t", quote = "", fill = TRUE)
#1.2读取KO数据----
eggnog.KEGG = read.table("E:/PCS/result/eggnog/eggnog.KEGG_ko.raw.txt",header = TRUE, sep = "\t", quote = "", fill = TRUE)
colnames(eggnog.KEGG) <- gsub("_clean", "", colnames(eggnog.KEGG))#删除原列名中"_clean"字符
#1.3ko与注释合并----
contig.kegg = merge(ko1_4,eggnog.KEGG,by.x = "KO",by.y = "KEGG_ko",all.y = F)#组装后kegg有重复，一个ko号对应多个代谢通路
length(unique(contig.kegg$PathwayL2))#查看PathwayL2有多少个单独项目
#1.4读取分类表----
sam = read.table("sampletab.txt",header = T)
#1.5读取SI----
SI.MID.MEAN.long=read.table("SI.MID.MEAN.LONG.txt",header = TRUE, sep = "\t", quote = "", fill = TRUE)


#2.数据清洗

gene_abun1=eggnog.KEGG

gene_abun1=as_tibble(gene_abun1)
gene_abun1=gene_abun1%>%
  select(-PS_1B,-PS_4B,-PS_8B)%>%
  column_to_rownames(var = "KEGG_ko") %>% # 将KEGG_ko设置为行名
  as.matrix()%>%                        # 将数据框转换为矩阵
  t()
gene_abun1=as_tibble(gene_abun1)
#head(eggnog.KEGG)[1:5,1:5]

#3.相关性分析
# 提取特定列1
first_columns <- gene_abun1

# 提取特定列2
last_column <- SI.MID.MEAN.long[,3] #均值

# 计算特定列1与特定列2的相关性
cor_results <- sapply(first_columns, function(x) {
  test <- cor.test(x, last_column)
  c(correlation = test$estimate, p_value = test$p.value)
})
cor_results_df <- as.data.frame(t(cor_results))# 打印相关性结果
#取显著的部分
cor_results_df.p=filter(cor_results_df,p_value<0.05)
cor_results_df.p$KO=rownames(cor_results_df.p)
#与分类、丰度表合并
targen=merge(cor_results_df.p,contig.kegg,by.x = "KO",by.y = "KO")
targen=targen%>%
  select(-PS_1B,-PS_4B,-PS_8B)
#4.基因处理
#删除人类疾病和无效数据
targen=filter(targen,PathwayL1!="Not Included in Pathway or Brite"&PathwayL1!="Human Diseases")
unique(targen$Pathway)
#转成长数据
targen_melt=melt(targen,id.vars =c("KO","correlation.cor","p_value","PathwayL1","PathwayL2","Pathway","KoDescription"))
#添加sam
targen_melt_sam=merge(targen_melt,sam,by.x = "variable",by.y = "sample")

ggplot(targen_melt_sam, aes(x = PathwayL1, y = value,fill = soiltreatmentabc)) +
  geom_bar(stat = "identity") +
  labs(x = "PathwayL2", y = "Value", title = "Bar Plot of PathwayL2 vs Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 如果 x 轴标签过长，可以旋转标签
ggplot(targen_melt_sam, aes(x = PathwayL1, y = value, fill = soiltreatmentabc)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "PathwayL2", y = "Value", title = "Bar Plot of PathwayL2 vs Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 创建柱状图并添加显著性标记
ggplot(targen_melt_sam, aes(x = PathwayL2, y = value, fill = soiltreatmentabc)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "PathwayL2", y = "Value", title = "Bar Plot of PathwayL2 vs Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(comparisons = list(c("healthC", "pathopoiesisC"), c("pathopoiesisC", "sterilizationC"), c("healthC", "sterilizationC")), 
              map_signif_level = TRUE, 
              position = position_dodge(0.9))
#unique(targen_melt_sam$soiltreatmentabc)
sum(is.na(targen_melt_sam))
targen_melt_sam$soiltreatmentabc <- as.factor(targen_melt_sam$soiltreatmentabc)

# 测试简单的 ggplot 和 ggsignif
ggplot(targen_melt_sam, aes(x = PathwayL1, y = value, fill = soiltreatmentabc)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  labs(x = "PathwayL2", y = "Value", title = "Bar Plot of PathwayL2 vs Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(comparisons = list(c("healthC", "pathopoiesisC")), 
              map_signif_level = TRUE, 
              position = position_dodge(0.9))
_
# 检查 PathwayL1 列和 value 列是否有缺失值
sum(is.na(targen_melt_sam$PathwayL1))
sum(is.na(targen_melt_sam$value))

# 简化数据集进行测试
targen_melt_sam_sample <- targen_melt_sam[1:10, ]

# 创建简单的 ggplot 并添加显著性标记
ggplot(targen_melt_sam_sample, aes(x = PathwayL1, y = value, fill = soiltreatmentabc)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  labs(x = "PathwayL2", y = "Value", title = "Bar Plot of PathwayL2 vs Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(comparisons = list(c("healthC", "pathopoiesisC")), 
              map_signif_level = TRUE, 
              position = position_dodge(0.9))

targen_melt_sam$PathwayL1 <- as.factor(targen_melt_sam$PathwayL1)
targen_melt_sam$soiltreatmentabc <- as.factor(targen_melt_sam$soiltreatmentabc)
ggplot(targen_melt_sam, aes(x = PathwayL1, y = value, fill = soiltreatmentabc)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  labs(x = "PathwayL2", y = "Value", title = "Bar Plot of PathwayL2 vs Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_signif(comparisons = list(c("healthC", "pathopoiesisC")), 
              map_signif_level = TRUE, 
              position = position_dodge(0.9))
__
# 安装并加载 ggpubr 包
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
}
library(ggpubr)

# 使用 ggpubr 创建柱状图并添加显著性标记
ggbarplot(targen_melt_sam, x = "PathwayL2", y = "value", fill = "stage", 
          position = position_dodge(0.9), add = "mean_se") +
  labs(x = "PathwayL2", y = "Value", title = "Bar Plot of PathwayL2 vs Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = stage), comparisons = list(c("healthC", "pathopoiesisC")),
                     label = "p.signif", method = "t.test")
__
# 计算每对比较的 p 值
compare <- targen_melt_sam %>%
  group_by(PathwayL1) %>%
  summarise(p_value = t.test(value[soiltreatmentabc == "healthC"], value[soiltreatmentabc == "pathopoiesisC"])$p.value)
# 计算每对比较的 p 值,数据有一定的参考价值
compare <- targen_melt_sam %>%
  group_by(PathwayL2) %>%
  summarise(p_value = t.test(value[soiltreatmentabc == "healthC"], value[soiltreatmentabc == "pathopoiesisC"])$p.value)
compare <- targen_melt_sam %>%
  group_by(PathwayL2) %>%
  summarise(p_value = t.test(value[soiltreatmentabc == "healthC"], value[soiltreatmentabc == "sterilizationC"])$p.value)
compare <- targen_melt_sam %>%
  group_by(PathwayL2) %>%
  summarise(p_value = t.test(value[soiltreatmentabc == "pathopoiesisC"], value[soiltreatmentabc == "sterilizationC"])$p.value)

# 计算每对比较的 p 值
compare <- targen_melt_sam %>%
  group_by(PathwayL2) %>%
  summarise(p_value = t.test(value[soiltreatmentabc == "healthC"], value[soiltreatmentabc == "pathopoiesisC"])$p.value)
#5.基因丰度热图绘制
##5.1数据清洗
targen_melt_sam

# 根据 Pathway 对 value 求和
pathway_sum <- targen_melt_sam %>%
  group_by(Pathway) %>%
  summarize(total_value = sum(value, na.rm = TRUE))

# 合并 pathway_sum 和 targen_melt_sam，按 Pathway 进行匹配
targen_melt_sam_updated <- targen_melt_sam %>%
  left_join(pathway_sum, by = "Pathway") %>%
  select(-KO, -KoDescription, -value)  # 删除 KO, KoDescription, value 列

hdata=targen_melt_sam_updated%>%select(variable,Pathway,total_value)%>%reshape2::dcast(.,Pathway~variable,mean)
#将第一列作为行名，并删除原第一列，再转成矩阵
hdata=hdata%>%
  column_to_rownames(var = "Pathway") %>% # 将KEGG_ko设置为行名
  as.matrix()                        # 将数据框转换为矩阵
mat=hdata

hdata1 <- targen_melt_sam %>%
  select(Pathway, correlation.cor,PathwayL1) %>%
  distinct()#%>% #基于选择的列删除重复行
#arrange(desc(correlation.cor))

Heatmap(mat)

__________________________________________________
#write.csv(targen_melt_sam,file = "result/targen.csv")
hdata=targen_melt_sam%>%select(variable,KO,value)%>%reshape2::dcast(.,KO~variable,mean)
hdata
#将第一列作为行名，并删除原第一列，再转成矩阵
hdata=hdata%>%
  column_to_rownames(var = "KO") %>% # 将KEGG_ko设置为行名
  as.matrix()                        # 将数据框转换为矩阵
mat=hdata
mat=log10(hdata+1)

hdata1 <- targen_melt_sam %>%
  select(KO, correlation.cor,PathwayL1) %>%
  distinct()#%>% #基于选择的列删除重复行
#arrange(desc(correlation.cor))

head(mat)
head(targen_melt_sam)
##5.2热图绘制
# 定义列的顺序
column_order = c("H_1C", "H_4C", "H_8C", 
                 "CH_1C", "CH_4C", "CH_8C", 
                 "OH_1C", "OH_4C", "OH_8C",
                 "PS_1C", "PS_4C", "PS_8C", 
                 "CPS_1C", "CPS_4C", "CPS_8C", 
                 "OPS_1C", "OPS_4C", "OPS_8C",
                 "S_1C", "S_4C", "S_8C", 
                 "CS_1C", "CS_4C", "CS_8C", 
                 "OS_1C", "OS_4C", "OS_8C")
# 假设数据矩阵 mat 已经存在，并且列名与 column_order 中的一致
# 先手动调整数据矩阵的列顺序，使其与 column_order 一致
mat_ordered <- mat[, column_order]
# 创建列注释
col_annotation <- HeatmapAnnotation(
  Group = factor(rep(c("heal", "path", "ste"), each = 9), levels = c("heal", "path", "ste"))
)
# 使用 column_split 分裂列
column_split <- rep(c("heal", "path", "ste"), each = 9)
# 创建行注释向量，使用 "+" 表示正相关，"-" 表示负相关
row_annotation_vector <- ifelse(hdata1$correlation.cor > 0, "+", "-")
# 创建 PathwayL1 行注释向量
pathway_annotation_vector <- hdata1$PathwayL1

# 创建双重行注释对象
row_annotation <- rowAnnotation(Sign = row_annotation_vector, PathwayL1 = pathway_annotation_vector)

# 绘制热图并添加双重行注释
Heatmap(mat_ordered,
        top_annotation = col_annotation, # 添加列注释
        right_annotation = row_annotation, # 添加双重行注释
        column_split = column_split,  # 使用 column_split 参数进行分裂
        cluster_columns = FALSE,  # 关闭列聚类以保持顺序
        cluster_row_slices = TRUE,  # 聚类行切片
        cluster_rows = TRUE)  # 聚类行
