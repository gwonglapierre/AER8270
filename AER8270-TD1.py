import numpy as np
import matplotlib.pyplot as plt


def géométrie(Nb_panneau, rayon):
    theta=np.linspace(0,2*np.pi,Nb_panneau+1)+np.pi/Nb_panneau
    x1=rayon*np.cos(theta)
    x2=rayon*np.sin(theta)
    
    return x1,x2

def Champ_vit(X,Y,X_source,Y_source,sigma):
    X_diff=X-X_source
    Y_diff=Y-Y_source
    R_sq=X_diff**2+Y_diff**2
    u=sigma/(2*np.pi)*(X_diff/R_sq)
    v=sigma/(2*np.pi)*(Y_diff/R_sq)
    
    return u,v

def calc_pressure_distribution(u,v,vel_inf):
    vel_mag= np.sqrt(u**2+v**2)
    cp=1-(vel_mag/vel_inf)**2
    
    return cp

def fonc():
    Nb_panneau=8 #Nombre de panneau
    rayon=1      #Rayon
    vel_inf=1    #Vitesse à l'infini
    sigma=1      #Force de la source
    cp_distribution=np.zeros(Nb_panneau)
    x_cyl,y_cyl=géométrie(Nb_panneau, rayon)
    
    for i in range(Nb_panneau):
        X_source, Y_source= x_cyl[i],y_cyl[i]
        u,v= Champ_vit(x_cyl[i], y_cyl[i], X_source, Y_source, sigma)
        cp_distribution[i]=calc_pressure_distribution(u, v, vel_inf)
        
    
    
    plt.plot(np.degrees(np.linspace(0,2*np.pi,Nb_panneau)), cp_distribution)
    plt.show()
    return
fonc()







