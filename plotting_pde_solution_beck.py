#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:18:41 2022

@author: beckettpierce
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(20)
u = np.zeros(20)
def steadystate(x,k,D,r):
    for j in range(len(x)):
        u[j] = (k/np.sqrt(4*D*r))**(-1*np.abs(x[j])*np.sqrt(r/D))
    return u

uhat = steadystate(x=x,k=100,D=6.45,r=.1)
#k in ug/ML
#D in 10^-6 cm^2/s
fig, ax = plt.subplots()
ax.plot(x, uhat)
ax.set(xlabel='x (10^-3 cm)', ylabel='drug concentration (ug/mL)',
       title='Steady State Drug Diffusion w/Constant Source & Reaction')
ax.grid()
