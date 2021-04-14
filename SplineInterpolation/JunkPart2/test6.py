
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
import numpy
import matplotlib.pyplot as plt
input=open('C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt','r')



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

    return x, y, z

if __name__ == "__main__":

    x, y, z = initializationData()

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')


    vertices = [list(zip(x,y,z))]

    poly = Poly3DCollection(vertices, alpha=0.8)

    ax.add_collection3d(poly)

    ax.set_xlim(0,5)

    ax.set_ylim(0,5)

    ax.set_zlim(0,5) 
    plt.plot()