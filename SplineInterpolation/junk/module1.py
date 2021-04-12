import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
import numpy

##import sys

##import matplotlib
#import matplotlib.pyplot as plt
##from matplotlib.ticker import MaxNLocator
#from matplotlib import cm
##from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.colors import LinearSegmentedColormap


##import numpy as np

##from scipy import array, newaxis

input=open('C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt','r')

def findSpline(valOx,valOy,x,y,dict):
    tempX = x.index(valOx)
    tempY = y.index(valOy)
    q= float(dict.get((2,2)))
    qq= float(dict.get((1,6)))
    qqq= float(dict.get((tempX,tempY)))
    return dict[(tempX,tempY)]

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

    #print(a,b,c,d,nX,nY,sep='\t',end='\n')

    lines=input.readlines()

    #print(len(lines))

    for item in lines:
        x.append(float((item.split(';')[0])))
        y.append(float((item.split(';')[1])))
        z.append(float((item.split(';')[2])))
        f.append(float((item.split(';')[3])))

    #print(len(z))
    #print(z[:10])

    #от 0 до 80
    spl={}
    for i in range(nX+1):
        for j in range(nY+1):
            spl[i,j] = z[i*(nX+1)+j]

    exactSol={}
    for i in range(nX+1):
        for j in range(nY+1):
            exactSol[i,j] = f[i*(nX+1)+j]
    #print(spl)

    newx = numpy.arange (a, b, nX)
    newy = numpy.arange (c, d, nY)
    xgrid, ygrid = numpy.meshgrid(newx,newy)
    #qwert=xgrid[len(xgrid)-1][len(xgrid[len(xgrid)-1])-1]
    #print(qwert)

    #for i in range(nX+1):
    #    for j in range(nY+1):
    #        zgrid.append(spl[i,j])
    #zgrid = findSpline(xgrid[len(xgrid)-1][len(xgrid[len(xgrid)-1])-1],ygrid[len(ygrid)-1][len(ygrid[len(ygrid)-1])-1],x,y,spl)
    zgrid =z
    return xgrid, ygrid, zgrid

if __name__ == "__main__":

    x, y, z = initializationData()

    fig = pylab.figure()
    axes = Axes3D(fig)

    axes.plot_surface(x, y, numpy.array(), rstride=2, cstride=2,cmap = LinearSegmentedColormap.from_list ("red_blue", ['b', 'w', 'g'], 256))#color='green')# cmap = cm.viridis)

    pylab.show()