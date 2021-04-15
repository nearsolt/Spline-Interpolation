import matplotlib.pyplot as plt

input=open('..\\SplineInterpolation1.txt','r')

if __name__ == "__main__":

    x=[]
    f=[]

    through1DerivType1=[]
    through1DerivType2=[]
    through1DerivType3=[]
    through1DerivType4=[]

    through2DerivType1=[]
    through2DerivType2=[]
    through2DerivType3=[]
    through2DerivType4=[]

    lines=input.readlines()
    print('number of rows: %s'%(len(lines)))
    for item in lines:
        x.append(float((item.split(';')[0])))
        f.append(float((item.split(';')[1])))
        through1DerivType1.append(float((item.split(';')[2])))
        through1DerivType2.append(float((item.split(';')[3])))
        through1DerivType3.append(float((item.split(';')[4])))
        through1DerivType4.append(float((item.split(';')[5])))
        through2DerivType1.append(float((item.split(';')[6])))
        through2DerivType2.append(float((item.split(';')[7])))
        through2DerivType3.append(float((item.split(';')[8])))
        through2DerivType4.append(float((item.split(';')[9])))

    input.close()

    fig, ax = plt.subplots()
    ax.plot(x,f, label = 'exact solution', color = 'black')
    ax.plot(x, through1DerivType1, label = 'type 1', color = 'red')
    ax.plot(x, through1DerivType2, label = 'type 2', color = 'green')
    ax.plot(x, through1DerivType3, label = 'type 3', color = 'yellow')
    ax.plot(x, through1DerivType4, label = 'type 4', color = 'blue')
    ax.legend(title ='Exact solution and 1st deriv',ncol = 2, edgecolor = '#e0e9ed')
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(x,f, label = 'exact solution', color = 'black')
    ax.plot(x, through2DerivType1, label = 'type 1', color = 'red')
    ax.plot(x, through2DerivType2, label = 'type 2', color = 'green')
    ax.plot(x, through2DerivType3, label = 'type 3', color = 'yellow')
    ax.plot(x, through2DerivType4, label = 'type 4', color = 'blue')
 
    ax.legend(title ='Exact solution and 2nd deriv',ncol = 2, edgecolor = '#e0e9ed')
    plt.show()

    
    fig, (ax1,ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
    plt.subplots_adjust()#hspace = 1
    fig.set_figheight(5)
    fig.set_figwidth(8)

    ax1.plot(x,f, 'black')
    ax1.plot(x, through1DerivType1, 'red')
    ax1.plot(x, through2DerivType1, 'blue')
    #ax1.set_title('Type 1')
    #ax1.set_xlabel('X-axis', fontweight ='bold')
    #ax1.set_ylabel('Y-axis', fontweight ='bold')

    ax2.plot(x,f, 'black')
    ax2.plot(x, through1DerivType2, 'red')
    ax2.plot(x, through2DerivType2, 'blue')
    #ax2.set_title('Type 2')
    #ax2.set_xlabel('X-axis', fontweight ='bold')
    #ax2.set_ylabel('Y-axis', fontweight ='bold')

    ax3.plot(x,f, 'black')
    ax3.plot(x, through1DerivType3, 'red')
    ax3.plot(x, through2DerivType3, 'blue')
    #ax3.set_title('Type 3')
    #ax3.set_xlabel('X-axis', fontweight ='bold')
    #ax3.set_ylabel('Y-axis', fontweight ='bold')

    ax4.plot(x,f, 'black')
    ax4.plot(x, through1DerivType4, 'red')
    ax4.plot(x, through2DerivType4, 'blue')
    #ax4.set_title('Type 4')
    #ax4.set_xlabel('X-axis', fontweight ='bold')
    #ax4.set_ylabel('Y-axis', fontweight ='bold')
    plt.show()