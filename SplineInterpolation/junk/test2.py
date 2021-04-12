
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
import numpy
input=open('C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt','r')

def initializationData():
    x=y=z=f=[]
    lines=input.readlines()
    print(len(lines))

    #for item in lines:
    #    x.append(float((item.split(';')[0])))
    #    y.append(float((item.split(';')[1])))
    #    z.append(float((item.split(';')[2])))
    #    f.append(float((item.split(';')[3])))

    #xgrid, ygrid = numpy.meshgrid(x, y)
    #zgrid = numpy.array([[z for i in range(9)] for j in range(9)])
    #print(zgrid.ndim)
    #print([x for i in range(9)])

    #return xgrid, ygrid, zgrid

initializationData()




#print('number of rows: %s'%(len(lines)))

#x, y, z = makeData()

#fig = pylab.figure()
#axes = Axes3D(fig)
##axes.plot(x,y,z)
#axes.plot_surface(x, y, z, rstride=2, cstride=2,color='green')#cmap = LinearSegmentedColormap.from_list ("red_blue", ['b', 'w', 'g'], 256))#color='green')# cmap = cm.RdBu)

#pylab.show()
