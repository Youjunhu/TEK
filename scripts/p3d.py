import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, pi
#from mayavi import mlab
import subprocess
import f90nml
import matplotlib.colors as mcolors

nml = dict(f90nml.read('input.nmlt'))
cn = dict(nml['control_nmlt'])
nrad = cn['nrad']


def my_plot(fn, sn):
    surf = np.loadtxt(fn)
    print('surf.shape=',surf.shape)

    m = nrad  # radial
    n =  surf[:,0].size//m # poloidal
    print('m,n=',m,n)

    r = surf[:,0].reshape(m,n)
    z = surf[:,1].reshape(m,n)
    field = surf[:,2].reshape(m,n)
    pfn =  surf[:,3].reshape(m,n)
    tfn =  surf[:,5].reshape(m,n)
    #add data at theta=+pi (which are identical to theta=-pi), so that no margine appears there:
    r = np.concatenate((r[:,:],r[:,0:1]),axis=1) 
    z = np.concatenate((z[:,:],z[:,0:1]),axis=1)
    field = np.concatenate((field[:,:],field[:,0:1]),axis=1)

    print('r.shape, z.shape, field.shape=', r.shape, z.shape, field.shape)
    print( 'field.min(),field.max(),=', field.min(), field.max())
    

    fig, ax = plt.subplots()
    jrad = m//2 - 89 # choos a radial location
    f1d = field[jrad,:]
    print('Radial location: rhot =', tfn[jrad,0]**0.5)

    #theta = [-pi + 2*pi/n*i for i in range(n+1)]
    R0 = 10
    theta = np.arctan2(z, r-R0)
    theta[:, -1] = pi
    theta = theta[jrad, :]
    print(theta.min(), theta.max())

    ax.plot(theta/pi, f1d, 'b.')
    ax.set_xlabel(r'$\theta/\pi$', fontsize=20)
    ax.text(0.05,0.9, sn, transform=ax.transAxes, fontsize=20)
    #plt.savefig(f'{fn}1d_theta.png',bbox_inches='tight')
    np.savetxt(fn+'_th.txt', np.c_[theta, f1d])
    plt.show()    

for it in [5,]:
    t=4000*it+1
    fn1 = f"poloidal_plane_t{t:06}Apara_nneq0"
    sn = r'$\delta A_{\parallel}$'
    my_plot(fn1, sn)

    fn2 = f"poloidal_plane_t{t:06}_nh001Apara"
    sn = r'$\delta A_{\parallel}$'
    my_plot(fn2, sn)

    fig, ax = plt.subplots()
    theta1 = np.loadtxt(fn1+'_th.txt')[:,0]
    f1 = np.loadtxt(fn1+'_th.txt')[:,1]
    theta2 = np.loadtxt(fn2+'_th.txt')[:,0]
    f2 = np.loadtxt(fn2+'_th.txt')[:,1]

    ax.plot(theta1/pi, f1, 'b.')
    ax.plot(theta2/pi, f2, '-k', lw=2)
    ax.plot(theta2/pi, -f2, '-k', lw=2)
    ax.set_xlabel(r'$\theta/\pi$', fontsize=20)
    ax.text(0.5,0.9, sn, transform=ax.transAxes, fontsize=20)
    plt.savefig(f'envolop_apara.png',bbox_inches='tight')
    plt.show()

    fn1 = f"poloidal_plane_t{t:06}Phi_nneq0"
    sn = r'$\delta \Phi$'
    my_plot(fn1, sn)

    
    fn2 = f"poloidal_plane_t{t:06}_nh001Phi"
    sn = r'$\delta \Phi$'
    my_plot(fn2, sn)

    fig, ax = plt.subplots()
    theta1 = np.loadtxt(fn1+'_th.txt')[:,0]
    f1 = np.loadtxt(fn1+'_th.txt')[:,1]
    theta2 = np.loadtxt(fn2+'_th.txt')[:,0]
    f2 = np.loadtxt(fn2+'_th.txt')[:,1]

    ax.plot(theta1/pi, f1, 'b.')
    ax.plot(theta2/pi, f2, '-k', lw=2)
    ax.plot(theta2/pi, -f2, '-k', lw=2)
    ax.set_xlabel(r'$\theta/\pi$', fontsize=20)
    ax.text(0.05,0.9, sn, transform=ax.transAxes, fontsize=20)
    plt.savefig(f'envolop_phi.png',bbox_inches='tight')
    plt.show()
