import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
def my_plot(fn, title):
    fig, ax = plt.subplots()


    x = np.loadtxt('profiles.txt')[:,1] # normalized toroidal flux
    x = x[0::4]
    x = x**0.5

    data= np.loadtxt(fn)
    t = data[:,0] # second
    t = t*1000 #to ms
    flux = data[:, 1:]

    t, x = np.meshgrid(t, x, indexing='ij')

    #c = ax.pcolormesh(t, x, flux, shading='auto')
    norm = mcolors.CenteredNorm()
    #c = ax.pcolormesh(t, x, flux, cmap='bwr', norm=norm, shading='auto')
    c = ax.pcolormesh(t, x, flux, cmap='seismic', norm=norm, shading='auto')
    cbar = fig.colorbar(c, ax=ax)
    cbar.ax.tick_params(labelsize=20)
    ax.set_title(title, fontsize=20)
    ax.set_xlabel(r'Time (ms)',fontsize=20)
    ax.set_ylabel(r'$\rho_t$',fontsize=20)
    ax.tick_params(axis='both', which='both', labelsize=15)
    plt.savefig(fn+'.png',bbox_inches='tight')


my_plot('phi.txt', title=r"$\delta \Phi$")
my_plot('apara.txt', title=r"$\delta A_{\parallel}$")
plt.show()
