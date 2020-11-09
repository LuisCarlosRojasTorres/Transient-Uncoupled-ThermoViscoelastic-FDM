# -*- coding: utf-8 -*-
"""
File: CS1fdm.py
Created on Wed Feb 12 14:35:01 2020
TRANSIENT STATE
Resolve através do Metodo das Diferencas Finitas o problema
termo-vicoelastico de geracao de temperatura numa laje submetida a 
carregamento cizalhante.
author: Luis Carlos A. Rojas Torres
email: luiscarlos.bsf@oceanica.ufrj.br
"""
import numpy as np
import FDM1D as fdm
from LinearViscoelastic import LinearViscoelastic
import output as out
import matplotlib.pyplot as plt

#EDITAR totalHoras= 20,40,60hr

#Problem DATA
TotalTime=20     #en horas
f=0.25          #frequency [Hz]
tau0=150000      #Shear Stress [Pa]
lamb=0.214       #Conductivity []
dc=1150*2200     #Density x Specific Heat
#BC
dT=0
T0=23

elastic_modulus = [38,20.14,14.3,9.97,7.25,5.15,3.4,2.1,1.8] #in MPa
relaxation_times =[0.58,3.13,18.72,125.41,1042,10942,159569,5215397]    #in seconds

#Viscoelastic object initialization
PU=LinearViscoelastic(elastic_modulus,relaxation_times)
#mesh configuration
lenght=0.1
numelem = 20
he=lenght/numelem
#Transient Configuration

print("frequency = ",f, "Hz")
print("Total Time = ",TotalTime, "hr")
dtime = 60   #en segundos


C1 = lamb*dtime/(dc*he**2)
C2 = dtime/dc
C3 = C2*2*np.pi*f*(tau0**2)

TotalTime*=3600/dtime


L=fdm.initL(numelem,he)
T=fdm.initT(numelem,T0)
print("   -Submitted")
T=fdm.explicitSol(L,PU,T,f,C1,C2,C3,dtime,TotalTime)
print("   -Completed")
plt.plot(L,T)
#out.print2arrays(L,T,'cs1_fdm_f05')
