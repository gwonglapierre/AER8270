import numpy as np 
from matplotlib import pyplot as plt
from TD1_fonctions import *


npan = 8

X,Y = cercle(1,npan)
x,y = disc(X,Y)
s = np.zeros(npan)
phi = np.zeros(npan)


for i in range(npan) :
    dx = X[i+1] - X[i]
    dy = Y[i+1] - Y[i]
    s[i] = np.sqrt(dx**2+dy**2)
    phi[i] = np.arctan2(dx,dy)


fig, ax1 = plt.subplots()
ax1.plot(X,Y,color='k')
ax1.scatter(x,y,12,'r')
ax1.scatter(X,Y,10,'k')
ax1.axis('equal')
ax1.grid()


plt.show()

