#! /usr/bin/env python3

import argparse
import gzip
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

def reads_count(file_text, dic_count):
    for i, line in enumerate(file_text) :
        if i % 4 == 1:
            if line.rstrip() in dic_count:
                dic_count[line.rstrip()] += 1
            else:
                dic_count[line.rstrip()] = 1
    return dic_count

def length_stat(seq_dict, filename):
    f = open(filename, 'w')
    # 定义两个列表分别用来存储序列长度和数量
    li_num = []
    li_len = []
    for index, seq in enumerate(seq_dict):
        num = seq_dict[seq]
        length = len(seq)
        li_num.append(num)
        li_len.append(length)
        f.write('>'+"seq"+str(index)+"_"+str(num)+"_"+str(length)+'\n')
        f.write(seq+'\n')
    f.close()
    return np.array([li_len, li_num])

# 传入array
def len_num_stat(len_num_array):
    len_min = min(len_num_array[0])
    len_max = max(len_num_array[0])
    unique_num = []
    total_num = []
    
    for i in range(len_min, len_max+1):
        n = 0
        s = 0
        for j, m in zip(len_num_array[0], len_num_array[1]):
            if j == i:
                n += 1
                s += m
        unique_num.append(n)
        total_num.append(s)
    
    return np.array([range(len_min, len_max+1), unique_num, total_num])

# plot histogram
def plot_bar(stat_array, filename):
    # 条形宽度
    bar_width = 0.3
    index_total = np.arange(len(stat_array[0][:13]))+2
    index_unique = index_total + bar_width
    plt.bar(index_total, height=stat_array[2][:13], width=bar_width, color="b", label="total_reads")
    plt.bar(index_unique, height=stat_array[1][:13], width=bar_width, color="g", label="unique_reads")
    plt.legend() # 显示图例
    plt.xticks(index_total + bar_width/2, stat_array[0][:13])
    plt.ylabel("reads count")
    plt.title("small RNA length num")
    plt.savefig(filename)


def store_dict(dic, out):
    f = open(out, 'w')
    json.dump(dic, f)
    f.close()


if __name__ == "__main__":
    # 
    #start = time.time()
    parser = argparse.ArgumentParser(description="将Clean data文件变成unique fasta文件")
    # add arguments
    parser.add_argument("-i", action="store", help="input file name")
    parser.add_argument("-o", action="store", help="output file prefix")
    #parser.add_argument("-f", action="store", help="input file format")
    args = parser.parse_args()
    
    file_text = gzip.open(args.i, 'rt')
    dic = {}
    # 使用字典统计small RNA的数目，以序列作key，数目作value
    seq_dict = reads_count(file_text, dic)

    # 创建文件夹
    curr_path = os.getcwd()
    target_dir = curr_path + "/unique_fasta/"
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
    
    path_file_pre = target_dir + args.o 
    file_array = length_stat(seq_dict, path_file_pre+".fasta")
    len_num_array = len_num_stat(file_array)
    np.savetxt(path_file_pre + ".txt",  np.transpose(len_num_array), fmt="%d", delimiter="\t")
    plot_bar(len_num_array, path_file_pre+".pdf")
