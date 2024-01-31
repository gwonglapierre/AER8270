import numpy as np 
from matplotlib import pyplot as plt 

def cercle(R, n) :
    theta = (np.pi*(1/n+1))-np.linspace(0,2*np.pi, n+1)
    X = R*np.cos(theta)
    Y = R*np.sin(theta)

    return X, Y

def disc(X,Y) :
    n = len(X)-1
    x = np.zeros(n)
    y = np.zeros(n)
    for i in range(n) :
        x[i] = (X[i]+X[i+1])/2
        y[i] = (Y[i]+Y[i+1])/2

    return x, y



