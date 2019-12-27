#! ~/miniconda3/bin/python3
import sys
 
if __name__ == "__main__":
  # 初始化一个names的字典，内容为空
  # 字典中为name和出现数量的键值对
  names = {}
  # sys.stdin是一个文件对象。 所有引用于file对象的方法，
  # 都可以应用于sys.stdin.
  for name in sys.stdin.readlines():
      # 每一行都有一个newline字符做结尾
      # 我们需要删除它
      name = name.strip()
      if name in names:
          names[name] += 1
      else:
          names[name] = 1
 
  # 迭代字典,
  # 输出名字，空格，接着是该名字出现的数量
  for name, count in names.iteritems():
      sys.stdout.write("%d\t%s\n" % (count, name))
