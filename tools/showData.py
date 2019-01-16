#!/usr/bin/python 
import os,sys,time
import numpy as np
import matplotlib.pyplot as plt
import re

def as_num(x):
    if x[0]=='N':
        print("NaN !!")
        return -1e10
    m=re.compile("E").split(x)
    return float(m[0])*10**int(m[1])

C_t = 360.515
# C_t = 37555.7641
A_1 = 1410435417.714647/(C_t**2)
A_2 = 0.007019709312750/C_t
qcon= 0.007019709312750

if __name__ == "__main__":
    if len(sys.argv)<2:
        print("one arguement is necessary!")
        sys.exit(0)

    filename = sys.argv[1]
    fin = open(filename,'r')
    lines = fin.readlines()
    print("lines:",len(lines),end='')

    data = []

    count = 0
    for line in lines:
        while line[0]==' ':
            line = line[1:]
        if count==0:
            data.append([str(x) for x in re.split(' |\n',line) if x!=''])
        else:
            data.append([as_num(x) for x in re.split(' ',line) if x!=''])
        count += 1

    print(data[0])
    print(" rows:",count)

    col_num = len(data[0])
    row_num = count

    values = [0 for i in range(col_num)]
    for col in range(col_num):
        msg = []
        for row in range(1,row_num):
            msg.append(data[row][col])
        values[col]=msg

    fig=plt.figure()
    for col in range(1,col_num):
        plt.subplot(col_num-1,1,col)
        #plt.figure(str(data[0][col]))
        plt.plot(values[0],values[col],'.-',ms=2)
        plt.xlabel(str(data[0][0]))
        plt.ylabel(str(data[0][col]))

    # plt.figure()
    # plt.subplot(3,1,1) # rho u^2 + P / r
    # ru2 = np.dot(np.dot(values[2],values[6]),values[6])
    # ru2p = ru2+np.array(values[1])/3
    # plt.plot(-np.diff(ru2p)/np.diff(values[0]))
    
    # plt.subplot(3,1,2) # 2 / r * rho u^2  
    # ru2x = 2*ru2/np.array(values[0])
    # plt.plot(-ru2x)

    # plt.subplot(3,1,3) # 
    # plt.plot(-np.diff(ru2p)/np.diff(values[0])-ru2x[1:])

    plt.show()

