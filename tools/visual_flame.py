from utils import *

def smooth(a,WSZ): 
    # a: NumPy 1-D array containing the data to be smoothed 
    # WSZ: smoothing window size needs, which must be odd number, 
    # as in the original MATLAB implementation 
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ  
    r = np.arange(1,WSZ-1,2) 
    start = np.cumsum(a[:WSZ-1])[::2]/r 
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1] 
    return np.concatenate(( start , out0, stop )) 

if __name__ == "__main__":
    if len(sys.argv)<2:
        print("one arguement is necessary!")
        sys.exit(0)
    
    start, end = 0, -1
    if len(sys.argv)>2:
        start = int(sys.argv[2])
    if len(sys.argv)>3:
        end = int(sys.argv[3])

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
    print(" rows:",count)
    print(data[0])

    col_num = len(data[0])
    row_num = count

    values = [0 for i in range(col_num)]
    for col in range(col_num):
        msg = []
        for row in range(1,row_num):
            msg.append(data[row][col])
        values[col]=msg[start:end]

    fig, axs = plt.subplots(col_num-1, 1, figsize=(4,4))
    for col in range(1,col_num):
        ax = axs[col-1]
        #plt.figure(str(data[0][col]))
        # if col==1:
        #     plt.plot(values[0],[int(val*100) for val in values[col]],'b.-',ms=2)
        # else:
        ax.plot(values[0],values[col],'b.-',ms=2)
        ax.set_xlabel(str(data[0][0]))
        ax.set_ylabel(str(data[0][col]))

    #values[1] = smooth(values[1], 5)
    #values[3] = smooth(values[3], 5)
    #values[2] = cdiff(values[1])/cdiff(values[0])
    #values[4] = values[2] - values[3]

    #for col in range(2,col_num):
        #plt.subplot(col_num-1,1,col)
        #plt.plot(values[0],values[col],'r.-',ms=2)
        #plt.xlabel(str(data[0][0]))
        #plt.ylabel(str(data[0][col]))

    plt.show()
