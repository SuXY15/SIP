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

if __name__ == "__main__":
    if len(sys.argv)<3:
        print("At least two arguement is necessary!")
        sys.exit(0)

    dt = float(sys.argv[1])
    fileList = sys.argv[2:]
    listNum = len(fileList)
    colsList = [2,4,5,7]
    colsNum = len(colsList)
    for i,filename in enumerate(fileList):
        time = float(filename[-8:-4])*dt
        fin = open(filename,'r')
        lines = fin.readlines()
        print(filename,"lines:",len(lines))

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

        col_num = len(data[0])
        row_num = count

        values = [0 for i in range(col_num)]
        for col in range(col_num):
            msg = []
            for row in range(1,row_num):
                msg.append(data[row][col])
            values[col]=msg

        for j,col in enumerate(colsList):
            print(i,j)
            plt.subplot(listNum,colsNum,j+1+i*colsNum)
            plt.plot(values[0],values[col],'.-',ms=2)

            plt.xlabel(str(data[0][0]))
            plt.ylabel(str(data[0][col]))
    plt.show()

