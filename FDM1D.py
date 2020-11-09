# -*- coding: utf-8 -*-
"""
File: FDM1D.py
Created on Sun Feb  9 11:15:40 2020

author: Luis Carlos A. Rojas Torres
email: luiscarlos.bsf@oceanica.ufrj.br
"""
import numpy as np
import matplotlib.pyplot as plt
from LinearViscoelastic import LinearViscoelastic

def initL(numelem,he):
    L=np.linspace(0,numelem*he,numelem+1)
    return L

def initT(numelem,y):
    T=np.zeros(numelem+1)
    for i in range(T.size):
        T[i]=y
    return T
    
def fx(PU,lamb,T,f,tau0):
    w=2*np.pi*f
    return 0.5*w*(tau0**2)*PU.get_Jpp(T,f)/lamb

    
def explicitSol(L,PU,T,f,C1,C2,C3,dtime,TotalTime):
    
    for i in range(0,int(TotalTime)):
        for j in range(1,T.size-1):
            T[j]=C1*(T[j+1]-2*T[j]+T[j-1])+C3*PU.get_Jpp(T[j],f)+T[j]
        dummy=T[1]
        T[0]=C1*(T[1]-2*T[0]+dummy)+C3*PU.get_Jpp(T[0],f)+T[0]
    
    #plt.plot(L,T)
    return T
    
