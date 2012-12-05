#!/usr/bin/python

f=open("center",'r')
cent=f.read()
f.close()
(A,B,C)=[float(t) for t in cent.split()]

f=open("COV",'r')
lines=f.readlines()
close=100000.0
chosen=-1
for one in lines:
    (x,y,z)=[float(t) for t in one[30:54].split()]
    dist=(A-x)*(A-x)+(B-y)*(B-y)+(C-z)*(C-z)
    if (dist<close):
        close=dist
        chosen=one

print chosen[30:54]

