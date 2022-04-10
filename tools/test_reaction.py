import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
from utils import *

def React(T9):
    T9A= T9/(1.+0.067*T9)
    P1 = T9A**(5./6.)/(T9**(3./2.))
    P2 = (np.exp(-0.010*T9A**4)+5.56e-3*np.exp(1.685*T9A**(2./3.)))
    Q = 1.26e27*P1*np.exp(-84.165/(T9A**(1./3.)))/P2
    print(T9A)
    print(P1)
    print(P2)
    return Q

def test_reaction():
    T9 = np.arange(0.1, 3.3, 0.1)
    Q = React(T9)
    Q2= 1.26e27*np.exp(-84.165/(T9**(1./3.)))
    plt.plot(T9, Q)
    plt.plot(T9, Q2)
    plt.show()

def test_qcon():
    Q = 13.931*1.6022e-13*6.022e23/24*1000*0.35
    E1 = EOS_E(3.5e10, 0.1e9) + Q
    T = np.arange(1.0, 7.5, 0.1)*1e9
    E2 = EOS_E(3.5e10, T)
    print("%12.5e"%E1)
    plt.plot(T, E2)
    plt.show()
    
def test_val():
    Tu = np.arange(0.1, 3.2, 0.1)*1e9
    Tb = Tu2Tb(3.5e10, Tu)
    theta = Tb/Tu
    Ze = 1/3.*(theta-1)/theta*pow(Ea/(Tb/1e9),1./3.)
    plt.plot(Tu, Tb/1e9); plt.plot(Tu, Ze)
    plt.xscale('log'); plt.xlabel(r'$T_{u}$')
    plt.legend([r'$T_{b9}$','Ze'])
    plt.show()
    print(Tu, Tb)

if __name__ == '__main__':
    # test_reaction()
    # test_qcon()
    test_val()