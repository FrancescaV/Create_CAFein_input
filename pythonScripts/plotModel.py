import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import math
import linecache
import re
from math import sqrt
from math import sin
from math import cos
from math import pow
from math import pi
from math import log
from math import log10
from math import fabs
import os
import glob


# Constants
G_cgs      =  6.67259e-8
Rsun_toCm  = 6.95660e10
Msun_toG   = 1.98855e+33
MjPerMsun    = 9.54265748e-4
Mj_toG     = 1.89813e30
Rj_toCm    = 7.1492e9
HubbleTime_yr = 13.8e9
secPerMin  = 60.0
minPerHours = 60.0
hoursPerDay = 24.0
Zsun = 0.02
secPerDay = secPerMin * minPerHours * hoursPerDay

minsPerYears = 60.0*24.0*365.242199
secPerYears = secPerMin*minsPerYears
plt.rcParams.update({'font.size': 8})

historyFile = "/Users/francesca/Northwestern/Research/create_models_for_tides_and_oscillations/createWDmodelsForOscillationCode_FromMESA_fixedUnitsIn_04022013/output/WD/0p23Msun/log1673/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_09222013.dat"
PSFILE = "/Users/francesca/Northwestern/Research/create_models_for_tides_and_oscillations/createWDmodelsForOscillationCode_FromMESA_fixedUnitsIn_04022013/output/WD/0p23Msun/log1673/plots/MESAtoCAFein_log1673_0p23Msun_Z0p02_WD_deltaX5eM5in_09222013.eps"

m        = np.loadtxt(historyFile, skiprows=4, usecols=(0,))
logrho   = np.loadtxt(historyFile, skiprows=4, usecols=(1,))
logP     = np.loadtxt(historyFile, skiprows=4, usecols=(2,))
logT     = np.loadtxt(historyFile, skiprows=4, usecols=(3,))
csi      = np.loadtxt(historyFile, skiprows=4, usecols=(4,))
Vg       = np.loadtxt(historyFile, skiprows=4, usecols=(5,))
Astar    = np.loadtxt(historyFile, skiprows=4, usecols=(6,))
U        = np.loadtxt(historyFile, skiprows=4, usecols=(7,))
c1       = np.loadtxt(historyFile, skiprows=4, usecols=(8,))
BVfreq   = np.loadtxt(historyFile, skiprows=4, usecols=(9,))
lamb_freq= np.loadtxt(historyFile, skiprows=4, usecols=(10,))
g        = np.loadtxt(historyFile, skiprows=4, usecols=(11,))
Lr       = np.loadtxt(historyFile, skiprows=4, usecols=(12,))
V        = np.loadtxt(historyFile, skiprows=4, usecols=(13,))
del_ad   = np.loadtxt(historyFile, skiprows=4, usecols=(14,))
delVar   = np.loadtxt(historyFile, skiprows=4, usecols=(15,))
Cp       = np.loadtxt(historyFile, skiprows=4, usecols=(16,))
v_t      = np.loadtxt(historyFile, skiprows=4, usecols=(17,))
kappa_s  = np.loadtxt(historyFile, skiprows=4, usecols=(18,))
eps_ad   = np.loadtxt(historyFile, skiprows=4, usecols=(19,))
eps_s    = np.loadtxt(historyFile, skiprows=4, usecols=(20,))
c2       = np.loadtxt(historyFile, skiprows=4, usecols=(21,))
c3       = np.loadtxt(historyFile, skiprows=4, usecols=(22,))
c4       = np.loadtxt(historyFile, skiprows=4, usecols=(23,))
dlnLr_dlnr        = np.loadtxt(historyFile, skiprows=4, usecols=(24,))
P_scale           = np.loadtxt(historyFile, skiprows=4, usecols=(25,))
StartRadSurfLayer = np.loadtxt(historyFile, skiprows=4, usecols=(26,))
EndRadSurfLayer   = np.loadtxt(historyFile, skiprows=4, usecols=(27,))
Gamma1            = np.loadtxt(historyFile, skiprows=4, usecols=(28,))
entropy           = np.loadtxt(historyFile, skiprows=4, usecols=(29,))
kappa             = np.loadtxt(historyFile, skiprows=4, usecols=(30,))
chiRho            = np.loadtxt(historyFile, skiprows=4, usecols=(31,))
chiT              = np.loadtxt(historyFile, skiprows=4, usecols=(32,))
UtoUsurf          = np.loadtxt(historyFile, skiprows=4, usecols=(33,))


plt.close('all')

fig, ((ax1, ax1R), (ax2, ax2R), (ax3, ax3R), (ax4, ax4R), (ax5, ax5R), (ax6, ax6R), (ax7, ax7R)) = plt.subplots(nrows=7, ncols=2, sharex='col', sharey='row')
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.95, wspace=0.0, hspace=0.1)

xmin1 = 0.0
xmax1 = 0.99
xmin2 = xmax1
xmax2 = 1.0

#################################
ax1.plot(csi, Vg, color='k')
ax1.set_xlim(xmin1, xmax1)

ax1R.plot(csi, Vg, color='k')
ax1R.set_xlim(xmin2, xmax2)

ax1.set_yscale('log')
ax1R.set_yscale('log')

ax1.set_ylabel(r"$V_g$")
#################################
ax2.plot(csi, Astar, color='k')
ax2.set_xlim(xmin1, xmax1)

ax2R.plot(csi, Astar, color='k')
ax2R.set_xlim(xmin2, xmax2)

ax2.set_yscale('log')
ax2R.set_yscale('log')

ax2.set_ylabel(r"$A^*$")
#################################
ax3.plot(csi, U, color='k')
ax3.set_xlim(xmin1, xmax1)

ax3R.plot(csi, U, color='k')
ax3R.set_xlim(xmin2, xmax2)

ax3.set_yscale('log')
ax3R.set_yscale('log')

ax3.set_ylabel(r"$U$")
#################################
ax4.plot(csi, c1, color='k')
ax4.set_xlim(xmin1, xmax1)

ax4R.plot(csi, c1, color='k')
ax4R.set_xlim(xmin2, xmax2)
ax4R.set_ylim(1e-2, 10)

ax4.set_yscale('log')
ax4R.set_yscale('log')

ax4.set_ylabel(r"$c_1$")

#################################
ax5.plot(csi, c2, color='k')
ax5.set_xlim(xmin1, xmax1)

ax5R.plot(csi, c2, color='k')
ax5R.set_xlim(xmin2, xmax2)

ax5.set_yscale('log')
ax5R.set_yscale('log')

ax5.set_ylabel(r"$c_2$")
#################################
ax6.plot(csi, c3, color='k')
ax6.set_xlim(xmin1, xmax1)

ax6R.plot(csi, c3, color='k')
ax6R.set_xlim(xmin2, xmax2)

ax6.set_yscale('log')
ax6R.set_yscale('log')

ax6.set_ylabel(r"$c_3$")


plt.savefig(PSFILE)
plt.close('all')
