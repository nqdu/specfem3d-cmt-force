import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d 

# load data
stnm = 'CE.K400.BXZ.semd'
d_force = np.loadtxt(f"OUTPUT_FILES.FORCE/{stnm}")
d_cmt = np.loadtxt(f"OUTPUT_FILES.CMT/{stnm}")
d_mix = np.loadtxt(f"OUTPUT_FILES/{stnm}")

# new time vector
tmin = max(d_force[1,0],d_cmt[1,0],d_mix[1,0])
tmax = min(d_force[-2,0],d_cmt[-2,0],d_mix[-2,0])
t = np.linspace(tmin,tmax,500)

# interpolate
s_force = interp1d(d_force[:,0],d_force[:,1])(t)
s_cmt = interp1d(d_cmt[:,0],d_cmt[:,1])(t)
s_mix = interp1d(d_mix[:,0],d_mix[:,1])(t)

# plot
plt.figure(1,figsize=(14,5))
plt.plot(t,s_force + s_cmt,label='cmt+force')
plt.plot(t,s_mix,label='mix',ls='--')
plt.legend()
plt.savefig("seismo.jpg")