import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import interpolate

pi   = np.pi
h_   = 6.6260693e-34
hbar = h_/pi/2.
c_   = 2.99792458e8
k_   = 1.3806488e-23
N_   = 500*6.0221413e23

rho_0 = 3.5e10
rho_b = 3.5e10
T_0 = 0.1e9
T_b = 3.2e9
S_c = 466
scale = 1.0
Ea  = 84.165**3

pi23= 3*pi*pi

def EOS_E(rho, T):
    pisq313 = (pi*pi*3.)**(1./3.)
    p1 = 3./4./rho*pisq313*hbar*c_*(rho*N_)**(4./3.)
    p2 = 1./2.*N_*pisq313**2/3/hbar/c_*(rho*N_)**(-1./3.)*(k_*T)**2.
    return p1+p2

E_0 = EOS_E(rho_0, 0.0)
E_T0= EOS_E(rho_0, T_0)
E_Tb= EOS_E(rho_b, T_b)

C_t = np.sqrt(E_T0)/S_c/scale
qcon= 0.904*0.35
A_1 = E_T0/C_t**2./S_c**2.
A_2 = qcon/C_t

def Tu2Tb(rho, Tu):
    Qc = 5.6e13 # J/kg
    Tadd = 9/2.*(pi23)**(-2/3)*Qc*k_**(-2)*hbar*c_*rho**(1/3)*N_**(-2/3)
    return np.sqrt(Tu**2 + Tadd)

def as_num(x):
    if x[0]=='N':
        print("NaN !!")
        return 0
    m=re.compile("E").split(x)
    return float(m[0])*10**int(m[1])

def cdiff(data):
    if len(data)==1:
        return np.array([0])
    return np.array([0]+list(np.diff(data)/2))+np.array(list(np.diff(data)/2)+[0])

def plot_contour(ax, r, v):
    x = np.linspace(-np.max(r), np.max(r), 256)
    y = np.linspace(-np.max(r), np.max(r), 256)
    Y, X = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    f = interpolate.interp1d(r, v, fill_value=np.nan, bounds_error=False)
    V = f(R.flatten()).reshape(R.shape)
    
    ax.contour(X, Y, V, cmap='Greys', vmin=-3, vmax=3)
    ax.axis("equal")
