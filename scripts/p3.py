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
    r=np.concatenate((r[:,:],r[:,0:1]),axis=1) 
    z=np.concatenate((z[:,:],z[:,0:1]),axis=1)
    field = np.concatenate((field[:,:],field[:,0:1]),axis=1)

    print('r.shape, z.shape, field.shape=', r.shape, z.shape, field.shape)
    print( 'field.min(),field.max(),=', field.min(), field.max())
    
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    c = ax.plot_surface(r, z, field,alpha=0.5)
    ax.set_xlabel('R(m)')
    ax.set_ylabel('Z(m)')
    ax.set_title(sn)
    plt.show()
    
    fig, ax = plt.subplots()
    #c = ax.pcolormesh(r, z, field, shading='gouraud',cmap="viridis")
    #c = ax.pcolormesh(r, z, field, shading='gouraud',cmap="gnuplot2_r", vmin=-0.0844, vmax=0.14)
    #c = ax.pcolormesh(r, z, field, shading='gouraud',cmap="gnuplot2_r")
    #norm = mcolors.CenteredNorm()
    #c = ax.pcolormesh(r, z, field, shading='auto',cmap="seismic", norm=norm)
    #c = ax.pcolormesh(r, z, field, shading='auto',cmap="seismic")
    c = ax.pcolormesh(r, z, field, shading='auto',cmap="YlGnBu")
    #c = ax.pcolormesh(r, z, field, shading='auto',cmap="PuBuGn")

    ax.plot(r[0,:], z[0,:], '--g' )
    ax.plot(r[-1,:], z[-1,:], '--g' )

    cbar = fig.colorbar(c, ax=ax)
    ax.text(0.02,0.92, sn, transform=ax.transAxes, fontsize=20)
    ax.set_xlabel('R(m)')
    ax.set_ylabel('Z(m)')
    #ax.text(0.8,0.9, 'TEK', transform=ax.transAxes, fontsize=15)
    #c = ax.pcolormesh(r, z, field, shading='nearest')\
    #c = ax.pcolormesh(field)
    #ax.text(0.98,0.95, f'n={ntor}',horizontalalignment='right', transform=ax.transAxes,fontsize=16)
    ax.set_aspect('equal', 'box')
    plt.savefig(f'{fn}.png',bbox_inches='tight')
    plt.show()

    fig, ax = plt.subplots()
    f1d = np.sum(field**2, axis=1)/field.shape[1]
    ax.plot(np.sqrt(tfn[:,0]), f1d)
    ax.set_xlabel(r'$\rho_t$', fontsize=20)
    ax.text(0.05,0.9, sn, transform=ax.transAxes, fontsize=20)
    plt.savefig(f'{fn}1d.png', bbox_inches='tight')
    plt.show()

    fig, ax = plt.subplots()
    ax.set_prop_cycle(color = ['c', 'm', 'y', 'k'],
                       ls = ["-", "--", "-.", ":"],
                       lw = [2, 2, 3, 4])
    harmonics = np.fft.fft(field, axis=1)/n
    for mh in range(9,13):
        ax.plot(np.sqrt(tfn[:,0]), 2*np.abs(harmonics[:,mh]), label=f'm={mh}')

    ax.set_xlabel(r'$\rho_t$', fontsize=20)
    ax.text(0.05,0.9, sn, transform=ax.transAxes, fontsize=20)
    ax.legend(fontsize=20)
    plt.savefig(f'{fn}1d_mharmonics.png', bbox_inches='tight')
    plt.show()

    

    fig, ax = plt.subplots()
    f1d = field[m//2,:]
    theta = [-pi + 2*pi/n*i for i in range(n+1)]
    ax.plot(theta, f1d, 'b.')
    ax.set_xlabel(r'$\theta$', fontsize=20)
    ax.text(0.05,0.9, sn, transform=ax.transAxes, fontsize=20)
    plt.savefig(f'{fn}1d_theta.png',bbox_inches='tight')
    np.savetxt(fn+'_th.txt', np.c_[theta, f1d])
    plt.show()    

for it in [5,]:
    t=4000*it+1

    fn = f"poloidal_plane_t{t:06}_nh001Apara"
    sn = r'$|\delta A_{\parallel n}|$'
    my_plot(fn, sn)
    
    fn = f"poloidal_plane_t{t:06}_nh001Phi"
    sn = r'$|\delta \Phi_n|$'
    my_plot(fn, sn)
