import numpy as np 
from matplotlib import pyplot as plt
from TD1_fonctions import *



X,Y = cercle(1,8)
x,y = disc(X,Y)
plt.plot(X,Y)
plt.scatter(x,y)
plt.show