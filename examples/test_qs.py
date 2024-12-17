import numpy as np
import astropy.units as u
from cycler import cycler

import matplotlib.pyplot as plt

import ASSolarPy.GXBox as gxb
import ASSolarPy.REOwrapper as REOwrapper
import ASSolarPy.AtmoProfile as atp
import ASSolarPy.Atmosphere as atm
import ASSolarPy.RATAN as rtn
import ASSolarPy.Dipole as dipole
import ASSolarPy.Coordinates as coord
import ASSolarPy.QuietSun as qs
from ASSolarPy.Show import Show
#from test_show import Show

import os
import sys

m = sys.modules[__name__]
this_path = os.path.dirname(m.__file__)

atmo = atp.AtmoProfile(model = atp.AtmoProfile.Avrett2008) # , detalization = 1e6
reo = REOwrapper.REOLibrary('s:/Projects/Physics687/ProgramD64/agsGeneralRadioEmission.dll')

color_cycle = cycler(color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

frequencies = np.hstack([[2.8], np.arange(4, 19, 2)])*1e9

RSun = 960

map_lim = 1.2
half_n = 1200
det_from = 1.0025
det_to = 1.004
det_n = 1000
map_step = 0.01
kf_show = 3

#map_lim = 1.2
#half_n = 1200
#det_from = 1.0025
#det_to = 1.004
#det_n = 1000
#map_step = 0.0002
#kf_show = 3

q = qs.QuietSun()
zirin = qs.QuietSun._QS_temperature(frequencies)

posx = np.linspace(0, map_lim, half_n)
flux = q.calc_freefree_cut(atmo, reo, posx, frequencies)
ax_cut = None
for kf, c in zip(range(flux.shape[0]), color_cycle):
    f = flux[kf, :] / flux[kf, 0] * zirin[kf]
    ax_cut = Show.plot('Half-Sun', posx, f, ax = ax_cut, color = c['color'], label = str(round(frequencies[kf]*1e-9, 1)) + ' GHz')
ax_cut.legend()
ax_cut.set_xlabel(r'Distance, $R_{Sun}$')
ax_cut.set_ylabel(r'Temperature, ${}^\circ K$')
ax_cut.set_title(r'Avrett 2008, temperature through disk center (half-Sun)')
I_center = flux[:, 0]

posx = np.linspace(det_from, det_to, det_n)
flux = q.calc_freefree_cut(atmo, reo, posx, frequencies)
ax_cut = None
for kf, c in zip(range(flux.shape[0]), color_cycle):
    f = flux[kf, :] / I_center[kf] * zirin[kf]
    ax_cut = Show.plot('Detailed', posx, f, ax = ax_cut, color = c['color'], label = str(round(frequencies[kf]*1e-9, 1)) + ' GHz')
ax_cut.legend()
ax_cut.set_xlabel(r'Distance, $R_{Sun}$')
ax_cut.set_ylabel(r'Temperature, ${}^\circ K$')
ax_cut.set_title(r'Avrett 2008, , temperature through disk center (detailed, distance/height around TR)')

secax = ax_cut.secondary_xaxis('top', functions=(lambda _: ((_-1)*u.solRad).to_value(u.Mm), lambda _: (_*u.Mm).to_value(u.solRad) + 1))
secax.set_xlabel(r'Height, $Mm$')
plt.show(block = False)

flux_map, pos_r = q.calc_freefree_map(atmo, reo, map_step, map_lim, frequencies)
kf = kf_show
map_kf = flux_map[kf,:,:]
map_kf = map_kf *zirin[kf] / map_kf[map_kf.shape[0]//2, map_kf.shape[1]//2]
im = Show.image('Flux @'+str(frequencies[kf]*1e-9)+  ' GHz', map_kf, cblabel = r'Temperature, ${}^\circ K$')
im.axes.set(xticks=np.linspace(0, map_kf.shape[0], 5), xticklabels=np.linspace(-1.3*960, 1.3*960, 5));
im.axes.set(yticks=np.linspace(0, map_kf.shape[1], 5), yticklabels=np.linspace(-1.3*960, 1.3*960, 5));
im.axes.set_xlabel(r'Distance, $arcsec$')
im.axes.set_ylabel(r'Distance, $arcsec$')
im.axes.set_title(r'Avrett 2008, radiomap @'+str(round(frequencies[kf]*1e-9, 1))+' GHz, $T_{center}$ = '+str(round(zirin[kf]))+' ${}^\circ K$')

flux_flat = q._map_flat(pos_r)

#xg, yg = np.meshgrid(pos_r, pos_r)
#ax = Show.surf('', xg, yg, flux_map[7,:,:])
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
## A StrMethodFormatter is used automatically
#ax.zaxis.set_major_formatter('{x:.02f}')

arc_step = np.array((map_step, map_step))*RSun
arc_base = - np.array((map_lim, map_lim))*RSun
scan_x = pos_r*RSun
ax_flat = None
for kf, c in zip(range(len(frequencies)), color_cycle):
    map_kf = np.squeeze(flux_map[kf,:,:]) * zirin[kf] / flux_map[kf, flux_map.shape[1]//2, flux_map.shape[2]//2]
    map_kf = reo._temp2fluxpixel(map_kf, frequencies[kf], arc_step[0])
    flat_kf = flux_flat * zirin[kf] / flux_flat[flux_flat.shape[0]//2, flux_flat.shape[1]//2]
    flat_kf = reo._temp2fluxpixel(flat_kf, frequencies[kf], arc_step[0])
    scan = rtn.RATAN._convolve(map_kf, frequencies[kf], arc_step, arc_base, mode = 3)
    ax_flat = Show.plot('Scans Comparison', scan_x, scan, ax = ax_flat, color = c['color'], label = str(round(frequencies[kf]*1e-9, 1)) + ' GHz')
    flat = rtn.RATAN._convolve(flat_kf, frequencies[kf], arc_step, arc_base, mode = 3)
    ax_flat = Show.plot('Scans Comparison', scan_x, flat, ax = ax_flat, color = c['color'], linestyle='dashed', linewidth = 1)
    print(sum(scan)*arc_step[0])
    pass
ax_flat.legend()
ax_flat.set_xlabel(r'Distance, $arcsec$')
ax_flat.set_ylabel(r'"Scan Intensity", $s.f.u./arcsec$')
ax_flat.set_title(r'Scans (solid = Avrett 2008, dashed = flat disk), $T_{center}$ = Zirin 1991')

plt.show()