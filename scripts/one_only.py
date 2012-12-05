#!/usr/bin/python

from sets import Set

f=open("repeat")
ori_list=[]
for line in f:
    ori_list+=[line.split()]
out=[]
for t in ori_list:
    t.sort()
    out+=['_'.join(t)]
out_list=list(Set(out))
out_list.sort()
for line in out_list:
    print line

