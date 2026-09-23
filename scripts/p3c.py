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

    

for it in [5,]:
    t=4000*it+1
    fn = f"poloidal_plane_t{t:06}Apara_neq0"
    sn = r'$|\delta A_{\parallel m}|$'
    my_plot(fn, sn)
    
    fn = f"poloidal_plane_t{t:06}Phi_neq0"
    sn = r'$|\delta \Phi_m|$'
    my_plot(fn, sn)

