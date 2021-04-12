import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
import numpy
import matplotlib.pyplot as plt

input=open('C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt','r')

def findSpline(X,Y,dict,x,y):
    tempX = x.index(X)
    tempY = y.index(Y)
    return dict.get[(tempX,tempY)]
    



def initializationData():
    x=[]
    y=[]
    z=[]
    f=[]

    initData=input.readline().split(';')
    a = float(initData[0])
    b = float(initData[1])
    c = float(initData[2])
    d = float(initData[3])
    nX = int(initData[4])
    nY = int(initData[5])

    lines=input.readlines()

    for item in lines:
        x.append(float((item.split(';')[0])))
        y.append(float((item.split(';')[1])))
        z.append(float((item.split(';')[2])))
        f.append(float((item.split(';')[3])))


    ##от 0 до 80
    spl={}
    for i in range(nX+1):
        for j in range(nY+1):
            spl[i,j] = z[i*(nX+1)+j]

    #exactSol={}
    #for i in range(nX+1):
    #    for j in range(nY+1):
    #        exactSol[i,j] = f[i*(nX+1)+j]
    ##print(spl)

    newx = numpy.linspace(a, b, nX)
    newy = numpy.linspace(c, d, nY)
    xgrid, ygrid = numpy.meshgrid(newx,newy)
    zgrid = findSpline(xgrid, ygrid,spl,x,y)

    return xgrid, ygrid, zgrid

if __name__ == "__main__":

    x, y, z = initializationData()
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x,y,z, rstride=1, cstride=1)
    #ax.plot_surface(numpy.array(x),numpy.array(y),numpy.array(z,0.0), rstride=1, cstride=1)
    plt.show()
    #fig = plt.figure()
    #axes = Axes3D(fig)

    #ax.plot_surface(x, y, z, rstride=2, cstride=2)#,cmap = LinearSegmentedColormap.from_list ("red_blue", ['b', 'w', 'g'], 256))#color='green')# cmap = cm.viridis)

    #plt.show()

