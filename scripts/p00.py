import matplotlib.pyplot as plt
import numpy as np

def my_plot(fn, sn, yscale='linear'):
    data = np.loadtxt(fn)
    time = data[:, 0]*1000 # from seconds to ms
    y = data[:, 1:]
    index = np.unravel_index(np.argmax(y, axis=None), y.shape)
    y_max = y[:,index[1]]
    print(index[1])
    tfn0 = np.loadtxt('profiles.txt')[index[1]*4,1] # normalized toroidal flux
    rhot0 = tfn0**0.5
    print(rhot0)
    np.savetxt(f'1d_{fn}', np.c_[time, y_max])

    fig, ax = plt.subplots()
    ax.plot(time, y_max, label='')
    ax.text(0.05, 0.9, sn+r' @$\rho_t$ ='+f'{rhot0:.2f}', transform=ax.transAxes, fontsize=20)
    #ax.legend(fontsize=20, loc='lower left')
    ax.set_yscale(yscale)
    new_xlim = [0, 0.75]

    ax.set_xlabel('Time (ms)', fontsize=20)
    
#    ax.text(0.05, 0.8, 'n=0', fontsize=20, transform=ax.transAxes)
    ax.tick_params(axis='both', which='both', labelsize=15)
    plt.savefig(f'{fn}{yscale}.png',bbox_inches='tight')



fn = 'phi.txt'
sn = r'$\delta \Phi$'
my_plot(fn, sn, 'linear')

fn = 'apara.txt'
sn = r'$\delta A_{\parallel}$'
my_plot(fn, sn, 'linear')
plt.show()
