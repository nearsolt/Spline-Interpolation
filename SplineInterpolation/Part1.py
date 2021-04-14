import matplotlib.pyplot as plt

input=open('C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation1.txt','r')

if __name__ == "__main__":
    x=[]
    f=[]

    s1=[]
    s2=[]
    s3=[]
    s4=[]
    s11=[]
    s12=[]
    s13=[]
    s14=[]

    lines=input.readlines()
    print('number of rows: %s'%(len(lines)))
    for item in lines:
        x.append(float((item.split(';')[0])))
        f.append(float((item.split(';')[1])))
        s1.append(float((item.split(';')[2])))
        s2.append(float((item.split(';')[3])))
        s3.append(float((item.split(';')[4])))
        s4.append(float((item.split(';')[5])))
        s11.append(float((item.split(';')[6])))
        s12.append(float((item.split(';')[7])))
        s13.append(float((item.split(';')[8])))
        s14.append(float((item.split(';')[9])))
    input.close()

    plt.plot(x,f, 'black')
    plt.plot(x, s1, 'red')
    plt.plot(x, s2, 'green')
    plt.plot(x, s3, 'yellow')
    plt.plot(x, s4, 'blue')
    plt.show()
    plt.plot(x,f, 'black')
    plt.plot(x, s11, 'red')
    plt.plot(x, s12, 'green')
    plt.plot(x, s13, 'yellow')
    plt.plot(x, s14, 'blue')
    plt.show()

    fig, (ax1,ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
    ax1.plot(x,f, 'black')
    ax1.plot(x, s1, 'red')
    ax1.plot(x, s11, 'blue')

    ax2.plot(x,f, 'black')
    ax2.plot(x, s2, 'red')
    ax2.plot(x, s12, 'blue')

    ax3.plot(x,f, 'black')
    ax3.plot(x, s3, 'red')
    ax3.plot(x, s13, 'blue')

    ax4.plot(x,f, 'black')
    ax4.plot(x, s4, 'red')
    ax4.plot(x, s14, 'blue')

    plt.show()

