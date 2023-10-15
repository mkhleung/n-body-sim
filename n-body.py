# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 21:44:10 2022

@author: mleung
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 20:48:14 2019

@author: michael
"""

# Third-party libraries
import numpy as np
from IPython.display import clear_output
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import copy 
import time
import pylab as pl
from IPython import display

#===================================== Define constants
AU = 149597870700
m_earth = 5.972*10**(24)

G = 6.674*(10**(-11))*m_earth/(AU**3)

n_body = 4

dt = 60*60*12 #- 1=sec, 60*60*24 = 1DAY

#m = np.random.randn(n_body)*m_earth
m = np.ones((n_body,1))*10000000#*1000000/n_body                         #*m_earth

dia = np.ones((n_body,1))*20


posi = np.random.randn(n_body, 3)*2

#turn on the following for the solar mass at center
posi[0,:] = [0,0,0]                         #set location of "sun"
m[0] = 1000000                        #set mass of "sun"
dia[0,0] = m[0]**(1/3)                #set displayed diameter of the sun



#posi[1,:] = [2,2,2]                         #set location of "earth""
vi = np.zeros((n_body,3))
F = np.zeros((n_body,n_body,3))

timerange = 500

max_orbit=10

pos_tensor = np.zeros((timerange,n_body,3))

for t in range(timerange):
    print(t)
    #plt.close()
    for i in range(1,n_body):     # for each body i calc the force
        for j in range(n_body): # calc force acting on i due to j.
            dx=posi[i,0]-posi[j,0]
            dy=posi[i,1]-posi[j,1]
            dz=posi[i,2]-posi[j,2]
            r = np.sqrt(dx**2+dy**2+dz**2)
            if i==j:
                F[i,j,:] = [0,0,0]
            else:
                Fr = G*m[j] / r**2 *(-np.sign(r))
                F[i,j,:] = [Fr*dx/r, Fr*dy/r, Fr*dz/r]
                    
    F_final = np.zeros((n_body,3))
    for k in range(3):              # loop thru each force
        for i in range(n_body):     # for each body i calc the sum of all forces
            F_final[i,k] = np.sum(F[i,:,k])
    
#    a_final = np.zeros((n_body,3))        
#    for i in range(n_body):   
#        a_final[i,:] = F_final[i,:]/m[i]
#    
    vi = vi + F_final*dt
    posi = posi + vi*dt
    
    for i in range(n_body):     # for each body i calc the sum of all forces
        for j in range(3):
            if np.abs(posi[i,j])>max_orbit:
                posi[i,j]=(max_orbit-0.5)*np.sign(posi[i,j])
                vi[i,j]=0
            pos_tensor[t,i,j] = posi[i,j]
                
    #pos_oAU = posi/AU
    #print(t)   
    #print(F[1,0,:])   
    #clear_output()
for t in range(timerange):
    bound = max_orbit
    plt.clf()    
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(-bound,bound)
    ax.set_ylim3d(-bound,bound)
    ax.set_zlim3d(-bound,bound)
    ax.set_facecolor((1.0, 1, 1))
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.title(t) 
    ax.scatter(pos_tensor[t,:,0], pos_tensor[t,:,1], pos_tensor[t,:,2], s=dia, c='b');
    plt.pause(0.01)
    plt.show()
    fig.show()
    

    