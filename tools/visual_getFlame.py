from utils import *

def diffMax(x, y):
    """ Find max position of dy/dx
    """
    pos = np.argmax()
    return pos, x[pos]

def curvMax(x, y, n=10, N=100, show=False):
    """ Find max position of dy/dx with high resolution
    """
    from scipy import interpolate
    l = len(x)-1
    dx, dy = cdiff(x), cdiff(y)
    pos = np.argmax(dy/dx)
    lpos = 0 if pos-n<0 else pos-n
    rpos = l if pos+n>l else pos+n
    x, dy = x[lpos:rpos], dy[lpos:rpos]
    x_new = np.linspace(x[0], x[-1], N)
    f = interpolate.interp1d(x, dy, 'quadratic')
    dy_new = f(x_new)
    if show:
        plt.plot(x, dy, 'ob')
        plt.plot(x_new, dy_new, 'r')
    return x_new[np.argmax(dy_new)]

if __name__ == "__main__":
    if len(sys.argv)<2:
        print("one arguement is necessary!")
        sys.exit(0)
    fileList = []
    path = sys.argv[1]
        
    for files in os.listdir(path):
        filename = os.path.join(path,files)
        if os.path.isfile(filename) and filename[-4:-1]=='.tx' and filename[-12:-8]=='data':
            fileList.append(filename)
    if len(fileList)<1:
        print("one data file is necessary!")
        sys.exit(0)

    fileList = sorted(fileList)[10:]
    flamePos = []
    for i,fileName in enumerate(fileList):
        fin = open(fileName,'r')
        lines = fin.readlines()
        X = []
        T = []
        count = 0
        for line in lines:
            beg = 0
            while line[beg]==' ':
                beg += 1
            line = line[beg:]
            if count>0:
                data = [as_num(x) for x in re.split(' ',line) if x!='']
                X.append(data[0])
                T.append(data[4])
            count += 1

        pos = curvMax(X, T, n=3, N=10000, show=False)
        flamePos.append(pos)
    
    plt.figure("Flame Position")
    plt.plot(flamePos)
    
    plt.figure("Flame velocity")
    plt.plot(cdiff(flamePos))
    
    plt.show()