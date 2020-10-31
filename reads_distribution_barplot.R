library(tidyverse)
# 自定义函数，绘制reads长度分布直方图
read_txt <- function(filename, main_title) {
  file1 <- read_csv(filename, col_names = c('sequence', 'num', 'len'))
  # 定义一个向量，用来储存长度字符串
  v_length <- vector()
  # 定义一个向量，用来储存每个长度的序列总数
  v_total_reads <- vector()
  # 定义一个向量，用来储存unique序列数目
  v_uniq_reads <- vector()
  # 统计18-30nt长度的reads数目
  for (i in c(18:30)) {
    v_length <- c(v_length, i) 
    v_total_reads <- c(v_total_reads, sum(file1 %>% filter(len == i) %>% select(num))) 
    v_uniq_reads  <- c(v_uniq_reads, nrow(file1 %>% filter(len == i) %>% select(num)))
  }
  p <- tibble(v_length, v_total_reads, v_uniq_reads) %>% 
    gather(key = 'type', value = 'nums', -v_length) %>% # 转换数据框，宽变长
    ggplot(aes(x=factor(v_length))) + geom_bar(aes(y =nums, fill=factor(type)), stat = "identity",position = "dodge") +
    xlab('Length(nt)') + ylab('Nums') + theme(axis.ticks = element_blank()) +
    labs(fill='Reads', title = main_title) + scale_fill_discrete(labels=c('total reads', 'unique reads')) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

# 文件名数组
list_name <- c(paste('trim_',paste('2a-','MOCK-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2a-','RGSV-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2a-','RRSV-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2ab-','MOCK-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2ab-','RGSV-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2ab-','RRSV-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2b-','MOCK-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2b-','RGSV-',2:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2b-','RGSV-','S1','_raw',sep = ''),sep = ''),
  paste('trim_',paste('2b-','RRSV-',1:2,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('2b-','RRSV-','S1','_raw',sep = ''),sep = ''),
  paste('trim_',paste('DJ-','MOCK-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('DJ-','RGSV-',1:3,'_raw',sep = ''),sep = ''),
  paste('trim_',paste('DJ-','RRSV-',1:3,'_raw',sep = ''),sep = ''))

# 遍历文件数组，批量绘制直方图
for (i in list_name) {
  file_name <- strsplit(i,'_')[[1]][2]
  p <- read_txt(paste(i,'.txt',sep = ''), file_name)
  ggsave(p,file=paste(file_name,".png",sep=""))
}
