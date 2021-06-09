import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import assoc_laguerre

def R_nl(r,Z,n,l):
    a0=1
    a=(2*Z/(n*a0))**((l+3/2))
    b=(math.factorial(n-l-1)/(2*n*math.factorial(n+l)))**(-1/2)
    #l^k_n
    # k=2l+1
    # n=n-l-1
    c=assoc_laguerre(r,n-l-1,2*l+1)
    d=(2*Z*r)/(n*a0)
    e=r**l * np.exp(-Z*r/(n*a0))
    return a*b*c*d*e

start=32
stop=35
r=np.linspace(start,stop)

Z=2
n=4
l=0
e_1=-.5


f1=1/r
f2=np.exp(-(-2*e_1*r))
f3=R_nl(r,Z,n,l)

plt.plot(r,f1)
plt.plot(r,f2)
plt.plot(r,f3)
label=['f1','f2','f3']
plt.legend(label)
plt.show()
