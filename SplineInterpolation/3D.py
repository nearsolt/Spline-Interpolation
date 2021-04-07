import matplotlib.pyplot as plt

input=open('C:\\Users\\1\\source\\repos\\Diplom\\SplineInterpolation2.txt','r')

x=[]
y=[]
z=[]
f=[]

lines=input.readlines()
for item in lines:
    x.append(float((item.split(';')[0])))
    y.append(float((item.split(';')[1])))
    z.append(float((item.split(';')[2])))
    f.append(float((item.split(';')[3])))

print('number of rows: %s'%(len(lines)))

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.plot(x, y, z, label='parametric curve')

fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, f, s=0.5,depthshade=True,c='red')
#plt.show()


fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, s=0.5,depthshade=True,c='green')
plt.show()
#fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True)
#ax1.plot(x,y,z, 'black')
#ax2.plot(x,y,f, 'red')
#plt.show()

input.close()