import numpy as np
import matplotlib.pyplot as plt

nx = 257
rho = np.linspace(0., 1.0, nx)

n0 = 1.44131e17
c0 = 0.49123
c1 = 0.298228
c2 = 0.198739
c3 = 0.521298
nf = n0*c3*np.exp(-c2/c1*np.tanh((rho-c0)/c2))
np.savetxt("nf.txt", np.c_[rho, nf])

for i in [100, 200, 300,  400, 500, 600, 700, 800]:
    tf = rho*0 + i
    np.savetxt(f"tf{i}.txt", np.c_[rho, tf])

ni = rho*0 + 2.0e19
#ne = rho*0 + 2.0e19
#ni = ne - nf
ti = rho*0 + 1.0
te = rho*0 + 1.0

ne = ni #+ nf # charge neutrality

np.savetxt("ni.txt", np.c_[rho, ni])
np.savetxt("ne1.txt", np.c_[rho, ne]) 
np.savetxt("ti.txt", np.c_[rho, ti])
np.savetxt("te.txt", np.c_[rho, te]) 

plt.plot(rho, -1/c1*(1-np.tanh((rho-c0)/c2)**2))
plt.show()
plt.plot(rho, nf)
plt.show()
