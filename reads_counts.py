import gzip
import sys
import re

def reads_count(file_text, dic_count):
    pattren = re.compile(r'[A|G|C|T|N]')
    for line in file_text:
        if pattren.match(line[0]):
            if line.rstrip() in dic_count:
                dic_count[line.rstrip()] += 1
            else:
                dic_count[line.rstrip()] = 1
    return dic_count


def length_stat(seq_dict, filename):
    f = open(filename, 'w')
    for seq, num in seq_dict.items():
        f.write(seq+','+str(num)+','+str(len(seq))+'\n')
    f.close()


if __name__ == "__main__":
    file_text = gzip.open(sys.argv[1], 'rt')
    dic = {}
    seq_dict = reads_count(file_text, dic)
    write_file_name = sys.argv[1].split('/')[-1].split('.')[0] + '.txt'
    length_stat(seq_dict, write_file_name)
