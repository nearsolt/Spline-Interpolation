#import pylab
#from mpl_toolkits.mplot3d import axes3d
#from matplotlib.colors import linearsegmentedcolormap
#from matplotlib import cm
#import numpy

#def makedata ():
#    x = numpy.arange (-10, 10, 0.1)
#    y = numpy.arange (-10, 10, 0.1)
#    xgrid, ygrid = numpy.meshgrid(x, y)

#    zgrid = numpy.sin (xgrid) * numpy.sin (ygrid) / (xgrid * ygrid)
#    return xgrid, ygrid, zgrid

#x, y, z = makedata()

#fig = pylab.figure()
#axes = axes3d(fig)

#axes.plot_surface(x, y, z, rstride=2, cstride=2,cmap = linearsegmentedcolormap.from_list ("red_blue", ['b', 'w', 'g'], 256))#color='green')# cmap = cm.viridis)

#pylab.show()

#My 3d graph

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np

figure = plt.figure()
axis = figure.add_subplot(111, projection = '3d')

x = [1,2,3,4,5,6,7,8,9,10]
y = [5,6,7,8,2,5,6,3,7,2]
z = np.array([[1,2,6,3,2,7,3,3,7,2],[1,2,6,3,2,7,3,3,7,2]])

axis.plot_wireframe(x, y, z)

axis.set_xlabel('x-axis')
axis.set_ylabel('y-axis')
axis.set_zlabel('z-axis')

plt.show()