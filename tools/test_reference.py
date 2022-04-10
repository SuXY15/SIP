import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
import re

pi   = np.pi              # Math PI
h_   = 6.62607015e-34     # Planck's constant  [m^2 kg s^-1] or [J s]
hbar = h_/pi/2.           # reduced Planck constant 
c_   = 2.99792458e8       # speed of light [m/s]
k_   = 1.3806488e-23      # Boltzmann constant [m^2 kg s^-2 K^-1]
NA   = 6.0221413e23       # Avogadro's number [mol^-1]
N_   = 500*NA             # Number density [#electron/kmol]x[mol^-1] ?

gamma = 4./3.             # adiabatic index
ρ0 = 3.5e10            # density of fresh state [kg/m^3]
ρb = 3.5e10            # density of burnt state [kg/m^3]
T0 = 0.1e9               # Temperature of fresh state [K]
Tb = 3.2e9               # Temperature of burnt state 
S_L = 466                 # Laminar flame thickness
Ea  = 84.165              # Activation energy [K^1/3]

pisq = pi*pi
def EoS_e(rho, T):
    pisq313 = (3.*pisq)**(1./3.)
    p1 = 3./4./rho*pisq313*hbar*c_*(rho*N_)**(4./3.)
    p2 = 1./2.*N_*pisq313**2/3/hbar/c_*(rho*N_)**(-1./3.)*(k_*T)**2.
    return p1 + p2, p1, p2

# check the EoS
def check_EoS():
    T_arr = np.linspace(0, Tb, 100)
    e_p1_arr = np.array([EoS_e(ρ0, T)[1] for T in T_arr])
    e_p2_arr = np.array([EoS_e(ρ0, T)[2] for T in T_arr])
    plt.figure(1)
    plt.plot(T_arr, e_p1_arr, '--', label='e_p1')
    plt.plot(T_arr, e_p2_arr, '--', label='e_p2')
    plt.plot(T_arr, e_p1_arr+e_p2_arr, '-', label='e_tot')
    plt.legend()
    plt.show()

e_0  = EoS_e(ρ0, 0.0)[0] # [J/kg] or [m^2/s^2]
e_T0 = EoS_e(ρ0, T0)[0]
e_Tb = EoS_e(ρb, Tb)[0]
p_T0 = ρ0*e_T0*(gamma-1)
p_Tb = ρ0*e_Tb*(gamma-1)
c_T0 = np.sqrt(gamma*p_T0/ρ0)
c_Tb = np.sqrt(gamma*p_Tb/ρ0)

u0 = np.sqrt(e_Tb)

Qhat = 13.931 * 1.6022*1e-13 * NA / 0.024 / u0**2
λhat = c_Tb / u0 * (e_Tb - e_T0) / u0**2
Dhat = c_Tb / u0
C_1  = 3./4.*(gamma-1)*(3*pisq)**(1./3.)*hbar*c_*N_**(4./3.)*ρ0**(1./3.)/u0**2
C_2  = 1./2.*(gamma-1)*(3*pisq)**(2./3.)/3./hbar/c_*N_**(2./3.)*ρ0**(-1./3.)/u0**2*k_**2*(Tb-T0)**2
Ehat = Ea*((Tb-T0)/1e9)**(-1./3.)
Ahat = c_Tb / u0 * np.exp( Ehat * (Tb/(Tb-T0))**(-1./3.) )

print("e_T0   = %.6e [m^2/s^2]" % e_T0)
print("p_T0   = %.6e [Pa]"      % p_T0)
print("c_T0   = %.6e [m/s]"     % c_T0)

print("e_Tb   = %.6e [m^2/s^2]" % e_Tb)
print("p_Tb   = %.6e [Pa]"      % p_Tb)
print("c_Tb   = %.6e [m/s]"     % c_Tb)

print("u0   = %.6e [m/s]"       % u0)
print("u0^2 = %.6e [m^2/s^2]"   % u0**2)

print("Qhat = %.6e"     % Qhat)
print("λhat = %.6e"     % λhat)
print("Dhat = %.6e / Le"% Dhat)
print("C1   = %.6e"     % C_1)
print("C2   = %.6e"     % C_2)
print("Ehat = %.6e"     % Ehat)
print("Ahat = %.6e"     % Ahat)