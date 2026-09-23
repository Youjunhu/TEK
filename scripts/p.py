import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lstsq
from scipy.signal import find_peaks

twopi = 2*np.pi

def spectrum(fn, sn):
    data = np.loadtxt(fn)
    print(data.shape)
    t =  data[:,0]
    s0 = data[:,1:]

    t_window = t[-1] - t[0]
    print('angular frequency domain resolution (kHz):', 1/t_window*twopi/1000)
    s = np.fft.fft(s0, axis=0)/s0.shape[0]
    ncut = 80 # set frequency cutoff
    s = s[:ncut,:]
    s = s.T

    print(s.shape)
    flat_idx = np.argmax(np.abs(s))
    index = (np.unravel_index(flat_idx, s.shape))
    print(fn, ': Angular frequency (kHz*rad)=', index[1]/t_window*twopi/1000)
    radial_index = index[0]
    print('Radial index of the maximum =', radial_index)
    omega = np.asarray([n*twopi/t_window for n in range(ncut)])
    omega = omega/1000
    x = np.sqrt(np.loadtxt('profiles.txt')[::4,1])

    X, Y = np.meshgrid(x, omega, indexing='ij')
    print(X.shape, Y.shape)

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    #im = ax.pcolormesh(X, Y, np.abs(s), shading='auto', cmap = 'viridis')
    im = ax.pcolormesh(X, Y, np.abs(s), shading='auto', cmap="YlGnBu")
    fig.colorbar(im, ax=ax)
    ax.set_title(sn, fontsize=20)
    ax.set_xlim(x.min(), x.max())

    data = np.loadtxt('continua.txt')
    x = data[:,2]
    y = data[:,5]*twopi
    y = np.ma.masked_where(y<1, y)
    ax.plot(x,y, '.', color='grey', label=r'Toroidal continua')
    ax.set_ylim(0,800)


    ax.set_xlabel(r"$\rho_t$", fontsize=20)
    ax.set_ylabel(r"$\omega$ (kHz $\cdot$ rad)", fontsize=20)
    #ax.text(0.62,35,"BAE gap", fontsize=20)
    ax.tick_params(axis='both', labelsize=20)
    #ax.legend(fontsize=15, loc=(0.65,0.5))
    plt.savefig(fn+".png", bbox_inches='tight')
    return radial_index


def growth_frequency(fn, radial_index):
  data = np.loadtxt(fn)
  ns = 200
  time = data[ns:,0]
  signal = data[ns:, radial_index]

  # Find peaks (local maxima)
  peaks, _ = find_peaks(signal, width=30)
  peak_times = time[peaks]
  peak_values = signal[peaks]
  print('number of peaks found =',peaks.size)
  if peaks.size>0:
     period = peak_times[-1] - peak_times[-2]
     frequency = 1/period
  else:
      frequency = np.nan
     
  # solve problem Ax = b using least square:
  A = np.empty((peaks.size, 2))
  for i, t in enumerate(peak_times):
    A[i,0] = t
    A[i,1] = 1
  b = np.log(peak_values)
  slope, c = lstsq(A, b)[0] # Slope is the growth rate

  print('gamma=', slope/1000, 'omega=',frequency*2*np.pi/1000)
  print('gamma/omega=', slope/1000/(frequency*2*np.pi/1000))

  fig, ax = plt.subplots()
  ax.plot(time*1000, np.log(np.clip(signal, 1e-12, None)))
  ax.plot(peak_times*1000, np.log(peak_values), 'ok')
  ax.plot(peak_times*1000, slope*peak_times + c)
  ax.set_xlabel('t[ms]',fontsize=20)
  plt.savefig(f'{fn}_growth.png')



fn = 'phi.txt'
radial_index = spectrum(fn, sn=r'$\delta\Phi$')
growth_frequency(fn, radial_index)


fn = 'apara.txt'
radial_index = spectrum(fn, sn=r'$\delta A_{\parallel}$')
growth_frequency(fn, radial_index)
plt.show()
