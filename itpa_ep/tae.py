import numpy as np
import matplotlib.pyplot as plt
q0= 1.71+ 0.16*0.5**2
R0 = 10
B = 3.0
mu0 = 4*np.pi*1.0e-7 #permeability in SI unit
mass_i = 1*1.660539066e-27 
rhom = 2.0e19*mass_i

vA = B/np.sqrt(mu0*rhom)
T_tae = 4*np.pi*q0*R0/vA
print('f_tae (kHz)=', 1/T_tae/1000)
print('w_tae (kHz)=', 2*np.pi*1/T_tae/1000)
print('dt (seconds) chosen  =', 0.05*T_tae)


keV =1.6022e-16   #unit J
Ti =1*keV
vti = np.sqrt(2*Ti/mass_i)
elementary_charge=1.6022e-19
charge = 1*elementary_charge
print('larmor radius (m)= ', mass_i*vti/(B*charge))
