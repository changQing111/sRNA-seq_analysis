library(tidyverse)
library(plyr)
library(Biostrings)

# 将序列反向互补
recv_comp <- function(sequence) {
  DNA.str <- DNAStringSet(sequence)
  temp <- reverseComplement(DNA.str)
  return(as.character(temp))
}
# 将比对文件的header分割,增加num列, length列
split_header <- function(align_file) {
  v_header <- unlist(str_split(align_file$header, "_"))  # 将header字符串分割并unlist为向量
  len <- length(v_header)
  align_file$num <- as.integer(v_header[c(1:len) %% 3 == 2]) 
  align_file$length <- as.integer(v_header[c(1:len) %% 3 == 0])  
  return(align_file)
}
# 添加正负链列
add_minus <- function(align_file) {
  minus <- vector('character',length = length(align_file$align_stat))
  minus[align_file$align_stat==0] <- "+"
  minus[align_file$align_stat==16] <- "-"
  align_file$minus <- minus
  return(align_file)
}

# 读入比对到病毒基因组上的文件
read_align_file <- function(filename, names=colname) {
  align_file <- read_tsv(filename, col_names = names)
  chain_logical = align_file$align_stat=="16"
  # 将比对到负链的序列反向互补
  align_file$sequence[chain_logical] <- recv_comp(align_file$sequence[chain_logical])
  align_file <- split_header(align_file)  # 增加num, length列
  align_file <- add_minus(align_file)     # 增加minus列
  align_file$num[align_file$minus=="-"] <- (-align_file$num[align_file$minus=="-"])
  # 添加5'碱基类型
  align_file$prefer <- str_sub(align_file$sequence[1:length(align_file$sequence)],1,1)
  return(align_file)
}
# 所有18-30nt的reads数目统计
all_reads_num_stat <- function(reads_num_len) {
  v_length <- vector(length = 13)
  v_length <- paste(18:30, "nt", sep = "")
  # 定义一个向量，用来储存每个长度的序列总数
  v_total_reads <- vector(length = 13)
  # 定义一个向量，用来储存unique序列数目
  #v_uniq_reads <- vector(length = 13)
  # 将长度18-30的序列提取出来
  for (i in c(18:30)) {
    v_total_reads[i-17] <- sum(reads_num_len %>% filter(length == i) %>% select(num))
   # v_uniq_reads[i-17]  <- nrow(reads_num_len %>% filter(len == i))
  }
  data1 <- data.frame(length=v_length, total_sRNAs=v_total_reads) 
  # , unique_sRNAs=v_uniq_reads
    #gather(key = 'type', value = 'nums', -length)
  return(data1)
}

# 长度统计
len_stat <- function(align_file) {
  v_num <-  vector("numeric", 13)
  v_a_ratio <- vector("double", 13)
  v_t_ratio <- vector("double", 13)
  v_c_ratio <- vector("double", 13)
  v_g_ratio <- vector("double", 13)
  for(i in 18:30) {
    all_num <- sum(abs(align_file %>% filter(length == i) %>%select(num)))
    v_num[i-17] <- all_num
    v_a_ratio[i-17] <- sum(abs(align_file %>% filter(length == i, prefer=="A")%>%select(num))) / all_num 
    v_t_ratio[i-17] <- sum(abs(align_file %>% filter(length == i, prefer=="T")%>%select(num))) / all_num
    v_c_ratio[i-17] <- sum(abs(align_file %>% filter(length == i, prefer=="C")%>%select(num))) / all_num
    v_g_ratio[i-17] <- sum(abs(align_file %>% filter(length == i, prefer=="G")%>%select(num))) / all_num
  }
  len_dist <- data.frame(length=paste(18:30, 'nt', sep = ''),
                         num=v_num,A=v_a_ratio,U=v_t_ratio,C=v_c_ratio, G=v_g_ratio)
  len_dist <- len_dist %>% gather(key = 'base', value = "percent", -c(length, num))
  len_dist$new_num <- len_dist$num*len_dist$percent
  return(len_dist)
}

# barplot
plot_hist <- function(align_file) {
  p <- align_file %>% ggplot(aes(x=length, y=new_num)) + geom_bar(aes(fill=base), stat = "identity", width = 0.75) 
    #scale_y_continuous(expand = c(0, 0), limits=c(0, 4000000), breaks = c(seq(0, 4000000, by=1000000))) 
  p <-  p + labs(x="",y="reads of vsiRNAs",  fill="") +  theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(1,1), legend.justification = c(1,1),
          axis.text.x = element_text(size = 10, colour = "black"),
          axis.title.y = element_text(size="14"))
  #ggsave(p, file=paste(main_title, "_barplot.pdf", sep = ""))
  return(p)
}
# 绘制vsiRNA在基因组上的分布
## 整理数据
size_num_stat <- function(align_file) {
  v_pos <- unique(align_file$pos) # all the unique position
  v_num <- vector(length = length(v_pos)*6)  # a position has eigth condition
  for (i in seq_along(v_pos)) {
    j <- i*6-5
    # 所有的正链数目
    all_positi <- align_file %>% filter(pos==v_pos[i], minus=="+") %>% select(num)
    if (nrow(all_positi)==0) {
      v_num[j] <- 0
    } else {v_num[j] <- sum(all_positi)}
    # 所有的负链数目
    all_nega <- align_file %>% filter(pos==v_pos[i], minus=="-") %>% select(num)
    if (nrow(all_nega)==0) {
      v_num[j+1] <- 0
    } else {v_num[j+1] <- sum(all_nega)}
    # 所有的+21nt
    pos_two_one <- align_file %>% filter(pos==v_pos[i], length==21, minus=="+") %>% select(num)
    if (nrow(pos_two_one)==0) {
      v_num[j+2] <- 0
    } else {v_num[j+2] <- sum(pos_two_one)}
    # 所有的-21nt
    nega_two_one <- align_file %>% filter(pos==v_pos[i], length==21, minus=="-") %>% select(num)
    if (nrow(nega_two_one)==0) {
      v_num[j+3] <- 0
    } else {v_num[j+3] <- sum(nega_two_one)}
    # 所有的+20nt
    pos_two_two <- align_file %>% filter(pos==v_pos[i], length==22, minus=="+") %>% select(num)
    if (nrow(pos_two_two)==0) {
      v_num[j+4] <- 0
    } else {v_num[j+4] <- sum(pos_two_two)}
    # 所有的-20nt
    neg_two_two <- align_file %>% filter(pos==v_pos[i], length==22, minus=="-") %>% select(num)
    if (nrow(neg_two_two)==0) {
      v_num[j+5] <- 0
    } else {v_num[j+5] <- sum(neg_two_two)}
  }
  
  v_size <- rep(c("+all", "-all", "+21nt", "-21nt", "+22nt","-22nt"), length(v_pos))
  pos <- vector(length=length(v_pos)*6)
  for (i in seq_along(v_pos)) {
    j <- i*6-5
    pos[j:(j+5)] <- rep(v_pos[i],6)
  }
  temp <- data.frame(position = pos, Nums=v_num, Size=v_size)
  temp$Size <- factor(temp$Size, levels = c("+all","+21nt","+22nt","-all","-21nt","-22nt"))
  return(temp)
}
