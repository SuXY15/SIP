#!/usr/bin/python 
from utils import *
    
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

    values = [0 for i in range(col_num)]
    for col in range(col_num):
        msg = []
        for row in range(1,row_num):
            msg.append(data[row][col])
        values[col]=msg

    fig, axs = plt.subplots(col_num-1, 1, figsize=(6,10))
    fig.subplots_adjust(bottom=0.05, top=0.95, left=0.15, right=0.95, hspace=0.3)
    for col in range(1,col_num):
        ax = axs[col-1]
        ax.plot(values[0],values[col],'.-',ms=2)
        ax.set_xlabel(str(data[0][0]))
        ax.set_ylabel(str(data[0][col]))
    
    fig, ax = plt.subplots(figsize=(6,6))
    plot_contour(ax, np.array(values[0]), cdiff(values[2])/cdiff(values[0]))
    
    plt.show()

