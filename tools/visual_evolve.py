#!/usr/bin/python 
from utils import *
from threading import Thread

global fignum

def run(fileList,jump=1):
    
    for i,filename in enumerate(fileList):
        if(i%jump==0) and filename.split('/')[-1][:4]=="data":
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
                if count == 0:
                    pass
                elif count == 1:
                    data.append([str(x).strip('\n').strip('"') for x in 
                        line.split(' ') if x!='' and x!="VARIABLES="])
                elif count == 2:
                    pass
                else:
                    data.append([as_num(x) for x in re.split(' ',line) if x!=''])
                count += 1
            print(" rows:", len(data[0]))
            col_num = len(data[0])
            row_num = len(data)
            
            if i==0:
                fig, axs = plt.subplots(col_num-1, 1, figsize=(6,10))
                fig.subplots_adjust(bottom=0.05, top=0.95, left=0.15, right=0.95, hspace=0.3)

            values = [0 for i in range(col_num)]
            for col in range(col_num):
                msg = []
                for row in range(1,row_num):
                    msg.append(data[row][col])
                values[col]=msg
            
            for col in range(1,col_num):
                ax = axs[col-1]
                ax.cla()
                ax.plot(values[0],values[col],'.-',ms=0.6,lw=0.4)
                ax.set_xlabel(str(data[0][0]))
                ax.set_ylabel(str(data[0][col]))
            plt.pause(0.001)
    plt.ioff()
    plt.show()

if __name__ == "__main__":
    fileList = []
    jump = 1
    path = sys.argv[1]
    
    if len(sys.argv)>2:
        jump = int(sys.argv[2])
        
    for files in os.listdir(path):
        filename = os.path.join(path,files)
        if os.path.isfile(filename) and filename[0:4]=='data':
            fileList.append(filename)
    if len(fileList)<1:
        print("one data file is necessary!")
        sys.exit(0)

    run(sorted(fileList), jump=jump)
