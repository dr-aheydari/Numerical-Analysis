#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 23:16:32 2019

@author: aliheydari
"""

import numpy as np
import time
from numpy import linalg as LA


import warnings
warnings.filterwarnings("ignore")


# set up the grid

h = 0.025; #also 0.1, 0.025

x_Mgrid = int(1 / h)
N = x_Mgrid;
y_Lgrid = int(1 / h)

print(x_Mgrid)
# the interior points
inter_x = x_Mgrid - 1
inter_y = y_Lgrid - 1
#a = 1.0e-2;

# Creating the mesh so we could find the boundary values
x = np.linspace(0, 1, x_Mgrid + 1  )
y = np.linspace(0, 1, y_Lgrid + 1 )

##### MATRICES ARE STORED ONLY FOR THE HEAT MAP TO TEST ACCURACY #####

X, Y = np.meshgrid(x, y)

##### MATRICES ARE STORED ONLY FOR THE HEAT MAP TO TEST ACCURACY #####
# initialization of
    # the temporary vector that holds the finite difference values of the grid
V = np.zeros_like(X)
# the actual solution just so that we could check our work
real_sol = np.zeros_like(X)
# the RHS of the poisson equation (just initialized here)
f = np.zeros_like(X)
U = np.zeros_like(X)
X0 = np.zeros_like(X)
XN = np.zeros_like(X)
Y0 = np.zeros_like(X)
YN = np.zeros_like(X)

#RHS = np.zeros_like(((int(1/h)*int(1/h)),1),dtype = float);

# Dirchlet boundary conditions : Actual solution is cos(x+2y)




RHS = np.zeros(((inter_x ) * (inter_x),1),dtype = float)


# given tolerance
tol=1.0e-7


# here we we vary values of omega to find the optimal value

#for q in range(0,5):
#    # to increment by 0.1

omega = 1.5 #+ q/10;
counter = 0
stop_res = 1.0

# calculating the normal derivatives, real solutions and F
for i in range(0,N + 1):
    for j in range(0,N + 1):
        
        f[i,j] = -13 * np.sin(3*X[i,j] + 2 * Y[i,j]);
        real_sol[i,j] = np.sin(3 * X[i,j] + 2 * Y[i,j])  
        X0[i,j] = -3 * np.cos(3*X[i,j] + 2 * Y[i,j])
        XN[i,j]= 3 * np.cos(3*X[i,j] + 2 * Y[i,j])
        Y0[i,j]=  2 * np.cos(3*X[i,j] + 2 * Y[i,j])
        YN[i,j]= -2 * np.cos(3*X[i,j] + 2 * Y[i,j])
        
          
        
        
t_start = time.clock()
   
for k in range(1,N):
    
    # line x = 0
    resid = 2 * U[0,k] - U[1,k] - 0.5  * U[0,k+1] + 0.5 * U[0,k-1] + \
    -1 * ( 0.5 * h**2 * f[0,k] + h * X0[0,k])
    U[0,k] = .5 * resid
    # the corner (0,0)
    U[0,0] = 0.5 * U[0,1] - 1/2 * U[1,0] + 0.25 * h**2 *f[0,0] + h * XN[0,0]
    
    # line x = N
    resid = 2 * U[N,k] - U[N,k] - 0.5  * U[N,k+1] + 0.5 * U[N,k-1] + -1 *( 0.5 * h**2*f[N,k] + \
    h * XN[N,k])
    U[N,k] = .5 * resid   
    
    # corner (N,N)
    U[N,N] = 0.5 * U[N,1] - 1/2 * U[1,0] + 0.25 * h**2 *f[N,N] + h * X0[N,N]   
    
    # line y = 0    
    resid = 2 * U[k,0] - U[k,1] - 0.5  * U[k+1,0] + 0.5 * U[k-1,0] + \
    -1 * ( 0.5 * h**2 * f[k,0] + h * YN[k,0])
    U[k,0] = .5 *  resid
    
    # Corner (0,N)
    U[0,N] = 0.5 * U[1,N] + 0.5 * U[0,N-1]+ 0.25 * h**2 * f[0,N] + h * YN[0,N];
    
    resid = 2 * U[k,N] - U[k,N] - 0.5  * U[k+1,N] + 0.5 * U[k-1,N] + \
    -1 *( 0.5 * h**2*f[k,N] + h * Y0[k,N])    
    U[k,N] = .5 * resid     
    # Corner (N,0)
    U[N,0] = 0.5 * U[N,1] + 0.5 * U[N-1,0]+ 0.25 * h**2 * f[N,0] + h * Y0[N,0];    

    
while (stop_res > tol):
    

    dVmax = 0.0
    bounter = 0;
    
    for i in range(1,inter_y + 1):
        
        for j in range(1, inter_x + 1):
            # real solution at each point 
                
            RHS[bounter] = h**2 * (f[i,j]);

            # residual 
            resid = (U[i,j-1] + U[i-1,j] + U[i,j+1] + U[i+1,j] - 4.0 * U[i,j]) \
            - h**2 * f[i,j]
            dV = 0.25 * omega * resid
            U[i,j]+= dV
            dVmax = np.max([np.abs(dV),dVmax])
#            print(bounter)
            bounter += 1;
            
    # calculating the total residual    
    stop_res = dVmax/np.max(np.abs(U))

    counter += 1
    
t_end = time.clock()

        

        
print ("SOR with Omega = {0} <----> CPU time = \t{1:0.2f} \t <----> iterations = {2}"\
       .format(omega,t_end - t_start,counter))
        
        
        
