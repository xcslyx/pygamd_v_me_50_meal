import numpy as np
import os
import re
from numba import jit


def EuclideanDistances(A, B):
    BT = B.transpose()
    vecProd = np.dot(A,BT)
    SqA =  A**2
    sumSqA = np.matrix(np.sum(SqA, axis=1))
    sumSqAEx = np.tile(sumSqA.transpose(), (1, vecProd.shape[1]))
    SqB = B**2
    sumSqB = np.sum(SqB, axis=1)
    sumSqBEx = np.tile(sumSqB, (vecProd.shape[0], 1))
    SqED = sumSqBEx + sumSqAEx - 2*vecProd
    SqED[SqED<0]=0.0
    ED = np.sqrt(SqED)
    return ED


@jit(nopython=True)
def jisuan(a,b):
    c=np.transpose(b)
    a+=b
    a+=c
    return a


si=open('size.txt','r')
size=si.read()
size=eval(size)

s2=[]
for i in range(len(size)):
    s1=[]
    for j in range(len(size)):
        s1.append(((size[i]+size[j])/2))    
    s2.append(s1)
sig=np.array(s2)



li=open('lambda.txt','r')
lize=li.read()
lize=eval(lize)

l2=[]
for i in range(len(lize)):
    l1=[]
    for j in range(len(lize)):
        l1.append((lize[i]+lize[j])/2)    
    l2.append(l1)
l=np.array(l2)

onepath = "/home/tsc/FUS/512/1weimiao/3/" 
threepath = "/home/tsc/FUS/512/1weimiao/energy/1/" 
zhipath ="/home/tsc/FUS/512/1weimiao/index_lianzhixin_r/"


total_txt = os.listdir(onepath)
num = len(total_txt)
list1 = range(num) 
files = os.listdir(onepath)
lj_eps=0.8368

aa=0
wenjian=1
ss=0
for a in list1: 
    name = total_txt[a]
    fp= open(onepath+name, 'r')
    xMat=fp.read()
    xMat=eval(xMat)
    mn=[]
    allcount=0
 
    for i in range(526):
        mn.append([])
        for j in range(526):
            mn[i].append(0)
    mn=np.array(mn)
    mn.dtype='float64'
    m= open(threepath+name, 'w')
    zhi=open(zhipath+name,'r')
    zhi=zhi.read()
    zz=eval(zhi)
    zz=sorted(zz,key=lambda x:x[2])
    print(zz[0]) 
    lian=np.arange(512)
    for jj in np.setdiff1d(lian,[zz[0][0],301]):
        c=EuclideanDistances(np.array(xMat[zz[0][0]]),np.array(xMat[jj]))
        c=np.array(c)
        if c.min()>4:
            ss+=1
        else:
            aa+=((4*lj_eps*((sig**12)/(c**12)-(sig**6)/(c**6))+lj_eps*(1-l))*(c<(2**(1/6)*sig)))
            aa+=((l*4*lj_eps*(((sig**12)/(c**12))-(sig**6)/(c**6)))*(c>=(2**(1/6)*sig))*(c<4))
            mn+=aa
            ss+=1
            aa=0
            print("第%d个第%d条链"%(wenjian,ss))
    wenjian+=1
    ss=0
    
    mn=np.array(mn)

    for i in range(526):
        for j in range(526):
                m.write(str(mn[i][j]))
                m.write(' ')
        m.write('\n')
    m.close()
     

