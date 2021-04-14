import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap

input=open('C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt','r')

if __name__ == "__main__":

    x=[]
    y=[]
    z=[]
    f=[]

    lines=input.readlines()
    print('number of rows: %s'%(len(lines)))
    for item in lines:
        x.append(float((item.split(';')[0])))
        y.append(float((item.split(';')[1])))
        z.append(float((item.split(';')[2])))
        f.append(float((item.split(';')[3])))
    input.close()
    fig1 = plt.figure(figsize=(9, 9))
    ax1 = fig1.add_subplot(111, projection='3d')
    surf1 = ax1.plot_trisurf(x, y, z, cmap = LinearSegmentedColormap.from_list('fstColor',['b', '#b37afa', '#6bfffd'],256))
    ax1.set_title('Spline')
    # Adding labels
    ax1.set_xlabel('X-axis', fontweight ='bold') 
    ax1.set_ylabel('Y-axis', fontweight ='bold') 
    ax1.set_zlabel('Z-axis', fontweight ='bold')

    fig2 = plt.figure(figsize=(9, 9))
    ax2 = fig2.add_subplot(111, projection='3d')
    surf2 = ax2.plot_trisurf(x, y, f, cmap = LinearSegmentedColormap.from_list('sndColor',['b', '#b37afa', '#6bfffd'],256))
    ax2.set_title('Exact solution')
    # Adding labels
    ax2.set_xlabel('X-axis', fontweight ='bold') 
    ax2.set_ylabel('Y-axis', fontweight ='bold') 
    ax2.set_zlabel('Z-axis', fontweight ='bold')

    plt.show() 