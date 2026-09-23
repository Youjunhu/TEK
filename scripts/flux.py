import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, pi
import subprocess

def plot(fn, sn):

    data= np.loadtxt(fn)
    nsub = data[0,:].size - 2 # radial gripoint number
    t1d = data[:,0]*1000 # time, from ms to second
    flux_av = data[:, 1] 
    flux = data[:, 2:2+nsub]
    if fn[0:4] == 'heat':
         flux = flux/10**3 #from W/m^2 to kW/m^2
         flux_av = flux_av/10**3
    else:
         flux = flux/10**17
         flux_av = flux_av/10**17
    
    coor = np.loadtxt('xgrid.txt')
    pfn = coor[:,0]
    tfn = coor[:,1]
    
    xmin = pfn[0]
    xmax = pfn[-1]
    dx = (xmax-xmin)/nsub
    x = np.linspace(xmin, xmax-dx, nsub)

    rhot = np.interp(x, pfn, tfn)**0.5
    t, rhot = np.meshgrid(t1d, rhot, indexing='ij')
    
    fig, ax = plt.subplots()
    #c = ax.pcolormesh(t, rhot, flux, shading='auto')
    c = ax.pcolormesh(t, rhot, flux, shading='auto', cmap="YlGnBu")
    cbar = fig.colorbar(c, ax=ax)
    cbar.ax.tick_params(labelsize=20)

    ax.set_title(sn, fontsize=20)
    ax.set_xlabel(r'Time (ms)',fontsize=20)
    ax.set_ylabel(r'$\rho_t$',fontsize=20)
    plt.savefig(fn+'.png',bbox_inches='tight')
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(t1d, flux_av, lw=3)
    ax.set_xlabel('Time (ms)', fontsize=20)
    ax.set_ylabel(sn, fontsize=20)
    #ax.text(0.1, 0.9, yscale+' scale', fontsize=20, transform=ax.transAxes)
    ax.tick_params(axis='both', which='both', labelsize=15)
    plt.savefig(f'{fn}_1d.png',bbox_inches='tight')
    plt.show()

    
fn = "heat_flux_ns1.txt"
plot(fn, sn = r"Heat flux ($kW/m^2$)")    

fn = "ptcl_flux_ns1.txt"
plot(fn, sn= r"Particle flux $(10^{17}s^{-1}/m^{2}$)")    


