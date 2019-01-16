#!/usr/bin/python 
import os,sys,time
import numpy as np
import matplotlib.pyplot as plt
import re
from threading import Thread

global fignum

def as_num(x):
    if x[0]=='a':
        print("NaN !!")
        return 0
    m=re.compile("E").split(x)
    return float(m[0])*10**int(m[1])

def run(fileList,jump=1):
    fig = plt.figure("displayData")
    ax = fig.add_subplot(111)

    for i,filename in enumerate(fileList):
        if(i%jump==0):
            fin = open(filename,'r')
            lines = fin.readlines()
            print(filename,"\tlines:",len(lines),end='')
            data = []
            count = 0
            for line in lines:
                beg = 0
                while line[beg]==' ':
                    beg += 1
                line = line[beg:]
                if count==0:
                    data.append([str(x) for x in line.split(' ') if x!=''])
                else:
                    data.append([as_num(x) for x in re.split(' ',line) if x!=''])
                count += 1
            print(" rows:",count)
            col_num = len(data[0])
            row_num = count
            values = [0 for i in range(col_num)]
            for col in range(col_num):
                msg = []
                for row in range(1,row_num):
                    msg.append(data[row][col])
                values[col]=msg
            plt.clf()
            for col in range(1,col_num):
                plt.subplot(col_num-1,1,col)
                plt.plot(values[0],values[col],'.-',ms=0.6,lw=0.4)
                plt.xlabel(str(data[0][0]))
                plt.ylabel(str(data[0][col]))
            plt.pause(0.001)


if __name__ == "__main__":
    fileList = []
    jump = 1
    if len(sys.argv)==2:

        jump = int(sys.argv[1])
    
    path = '.'
    for files in os.listdir(path):
        filename = os.path.join(path,files)
        if os.path.isfile(filename) and filename[-4:-1]=='.tx':
            fileList.append(filename)
    if len(fileList)<1:
        print("one data file is necessary!")
        sys.exit(0)

    run(sorted(fileList), jump=jump)
