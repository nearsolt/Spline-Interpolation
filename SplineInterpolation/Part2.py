import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap as LSC

input=open('..\\SplinesOfTwoVariables.txt','r')

if __name__ == "__main__":

    x=[]
    y=[]
    f=[]
    m1=[]
    M2=[]

    lines=input.readlines()
    print('number of rows: %s'%(len(lines)))
    for item in lines:
        x.append(float((item.split(';')[0])))
        y.append(float((item.split(';')[1])))
        f.append(float((item.split(';')[2])))
        m1.append(float((item.split(';')[3])))
        M2.append(float((item.split(';')[3])))
    input.close()
    
    fig1 = plt.figure(figsize=(7, 7))
    ax1 = fig1.add_subplot(111, projection='3d')
    surf1 = ax1.plot_trisurf(x, y, f, cmap = LSC.from_list('fstColor',['#31053e', '#0a0753','#1e48db', '#0d6b04', '#edde1b','#ff1616'],256))
    ax1.set_title('The exact solution')
    # Adding labels
    ax1.set_xlabel('X-axis', fontweight ='bold') 
    ax1.set_ylabel('Y-axis', fontweight ='bold') 
    ax1.set_zlabel('Z-axis', fontweight ='bold')
    #plt.show() 
     
    fig2 = plt.figure(figsize=(7, 7))
    ax2 = fig2.add_subplot(111, projection='3d')
    surf2 = ax2.plot_trisurf(x, y, m1, cmap = cm.gnuplot)
    ax2.set_title('The representation of the spline through the first derivative')
    # Adding labels
    ax2.set_xlabel('X-axis', fontweight ='bold') 
    ax2.set_ylabel('Y-axis', fontweight ='bold') 
    ax2.set_zlabel('Z-axis', fontweight ='bold')
    #plt.show() 
    
    fig3 = plt.figure(figsize=(7, 7))
    ax3 = fig3.add_subplot(111, projection='3d')
    surf3 = ax3.plot_trisurf(x, y, M2, cmap = cm.gnuplot2)
    ax3.set_title('The representation of the spline through the second derivative')
    # Adding labels
    ax3.set_xlabel('X-axis', fontweight ='bold') 
    ax3.set_ylabel('Y-axis', fontweight ='bold') 
    ax3.set_zlabel('Z-axis', fontweight ='bold')
    plt.show() 