#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd


direction='y'
x0='x_dir*_'+direction+'*res_0*'
x1='x_dir*_'+direction+'*res_1*'
x2='x_dir*_'+direction+'*res_2*'
x3='x_dir*_'+direction+'*res_3*'
x4='x_dir*_'+direction+'*res_4*'
x_states=[x0,x1,x2,x3,x4]
x_lists=[]
for x in x_states:
    x_lists.append(glob.glob(x))
    x_lists[-1].reverse()

sta_num=0
orb_num=0
orb_data={}
for state in x_lists:
    orb_num=0
    
    for orbital_file in state:
        name='x_'+str(sta_num)+'_'+str(orb_num)
        orb_num=orb_num+1
        orb_data[name]=np.genfromtxt(orbital_file,dtype=float)
    sta_num=sta_num+1

state=0
orbitals=range(5)
n_cols=5
n_rows=1
fig,axes=plt.subplots(n_rows,n_cols)
for row,o in zip(axes,orbitals):
    key='x_'+str(state)+'_'+str(o)
    print(key)
    r=orb_data[key][:,0]
    phi=orb_data[key][:,1]
    r_range=np.argwhere(r>25)
    r=r[r_range]
    phi=phi[r_range]
    row.plot(r,phi)
    plt.title(key)

