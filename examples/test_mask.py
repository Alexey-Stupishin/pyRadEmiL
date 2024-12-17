import numpy as np
import astropy.units as u

import ASSolarPy.GXBox as gxb
import ASSolarPy.REOwrapper as REOwrapper
import ASSolarPy.AtmoProfile as atp
import ASSolarPy.Atmosphere as atm
import ASSolarPy.RATAN as rtn
import ASSolarPy.Dipole as dipole
import ASSolarPy.Coordinates as coord
import ASSolarPy.QuietSun as qs
from test_show import Show

import os
import sys

m = sys.modules[__name__]
this_path = os.path.dirname(m.__file__)

x = [1,2,3]
y = [4,5,6]

xv, yv = np.meshgrid(x, y)

#flux00t = [[1,2,3], [4,5,6], [7,8,9], [10,11,12]]
#flux00 = np.array(flux00t)
#print(flux00[1,2])
#sp = flux00.shape
#fflip = np.flip(flux00, 0)
#spf = fflip.shape


#flux0 = np.vstack([fflip, flux00[1::, :]])


atmo = atp.AtmoProfile(model = atp.AtmoProfile.Avrett2008)
reo = REOwrapper.REOLibrary('s:/Projects/Physics687/ProgramD64/agsGeneralRadioEmission.dll')

q = qs.QuietSun()
#frequencies = np.arange(4, 18, 1)*1e9
#posx = np.linspace(1.0025, 1.004, 1000)
frequencies = np.arange(4, 20, 4)*1e9
posx = np.linspace(1.0025, 1.004, 5)
flux = q.calc_freefree_cut(atmo, reo, posx, frequencies)

pos_km = ((posx-1)*u.solRad).to_value(u.Mm)
ax_cut = None
for kf in range(flux.shape[0]):
    ax_cut = Show.plot('', pos_km, flux[kf, :], ax = ax_cut)

map_step = 0.01
map_lim = 1.2
flux_map = q.calc_freefree_map(atmo, reo, map_step, map_lim, frequencies)
Show.image('flux1', flux_map[0,:,:])
Show.image('flux2', flux_map[-1,:,:])

scan_x = np.arange(-map_lim, map_lim-map_step, map_step)*960
ax_scan = None
for kf in range(len(frequencies)):
    scan = rtn.RATAN._convolve(np.squeeze(flux_map[kf,:,:]), frequencies[kf], np.array((map_step, map_step))*960, -np.array((map_lim, map_lim))*960)
    ax_scan = Show.plot('scans', scan_x, scan, ax = ax_scan)

pass

#a = np.array((1, 2, -3, 4, 5, 6, -7, 8))
#aa = np.sqrt(a)
#v = ~np.isnan(aa)

#h, t, d = a.los(0, 0)
#h, t, d = a.los(0.5, 0.5)
#h, t, d = a.los(0, 1.001)
#h, t, d = a.los(0, 1.003)
#h, t, d = a.los(0, 1.01)
#h, t, d = a.los(0, 2)

#v = atp._barometric_coef(1e8, [2e8, 3e8, 4e8], [1e6, 1e6, 1e6])

#a = np.array((1, 2, -3, 4, 5, 6, -7, 8))
#f = a >= 0
#aa = a[a >= 0]
#v = np.sqrt(a)

#v = np.sqrt(-1)
#p = np.isnan(v)

#xmap = np.array([[1,2], [3,4], [5,6], [7,8]])
#xbeam = np.array([1, 2])
#a = xmap.dot(xbeam)

#a0 = np.array([1, 2, 3, 4, 5])
#aa = np.atleast_2d(a0).T
## aa = np.reshape(a, np.flip(a.shape[0]))
#b = np.repeat(aa, 3, axis = 1)

#cc = coord.Coordinates(xy = [0.3, 0.6])
#lonlat = cc.lonlat
#dircos = cc.dir_cos

#d = dipole.Dipole()

# load sample box
box = gxb.GXBox(this_path + '/../Data/11312_hmi.M_720s.20111010_085818.W120N23CR.CEA.NAS_ext.sav')

# create details mask for this box
details_mask = box.get_details_mask
Show.image('mask', details_mask.astype(np.float64))

# create atmosphere model
atmo = atm.Atmosphere()

# set atmosphere profiles
# set common profile (for all mask values except specified below)
pressure = 3e15
heights =     [1,   2e8,   2.3e8, 3.0e8, 5.0e9]
temperature = [1e4, 1e4,   1.0e6, 1.2e6, 2.5e6]
density = pressure/np.asarray(temperature)
atmo.atmosphere_init(heights, temperature, density)
# set profile for the penumbra (mask == 6)
heights =     [1,   1.5e8, 2.5e8, 5.0e9]
temperature = [1e4,   1e4, 1.0e6, 2.0e6]
density = pressure/np.asarray(temperature)
atmo.atmosphere_mask_set(6, heights, temperature, density)
# set profile for the umbra (mask == 7)
heights =     [1,   1.0e8, 1.3e8, 2.0e8, 5.0e9]
temperature = [1e4,   1e4, 0.7e6, 0.9e6, 2.0e6]
density = pressure/np.asarray(temperature)
atmo.atmosphere_mask_set(7, heights, temperature, density)
# create atmosphere model by mask and profiles
# mask_model = atmo.phys_atmosphere_mask(details_mask)

# load radioemission library
reo = REOwrapper.REOLibrary('s:/Projects/Physics687/ProgramD64/agsGeneralRadioEmission.dll')
# reo = REOwrapper.REOLibrary(this_path + '/../bin/agsGeneralRadioEmission.dll')
# show version
print(reo.get_version)

# prepare library to calculate box

# get RATAN position angle
ratan = rtn.RATAN(this_path + '/../Data/20111010_130044_sun0_out.fits', type = rtn.RATAN.FILE_FITS)
pos_angle = ratan.position_angle

beam = ratan.beam(5e9)

map_step = 1
freefree = 0
res = reo.prepare(box, map_step, pos_angle, freefree = 0)

# apply atmosphere model to calculating box
res_atm = reo.set_atmosphere(atmo, mask = details_mask)

# calculate radiomaps and scans
frequency = 5.7e9
res_calc = reo.map_calculate(frequency, mode = rtn.RATAN.BEAM_2023)

# plot some results
Show.image_log('Left Tau', res_calc['tauL'])
scanL = res_calc['scanL']
x_scan = np.linspace(res_calc['scan_info']['scan_lim'][0], res_calc['scan_info']['scan_lim'][1], num = res_calc['scan_info']['scan_size'])
ax_scan = Show.plot('scanL', x_scan, res_calc['scanL'])

scan = rtn.RATAN._convolve(res_calc['fluxL'], frequency, [map_step, map_step],  res['arc_box'][:, 0])
Show.plot('', x_scan, scan, ax = ax_scan)

Show.image_log('Left Flux', res_calc['fluxL'])
# convert fluxes to temperature
Show.image_log('Left Flux in K', reo.fluxpixel2temp(res_calc['fluxL'], frequency))

# expand scan limits and convolve again with another diagram
scan_lim = res['arc_box'][0, :] + [-40.5, 40.7]
scan_data = reo.map_convolve(res_calc['fluxL'], 5.7e9, mode = rtn.RATAN.BEAM_STOH, scan_lim = scan_lim);
x_scan = np.linspace(scan_data['scan_info']['scan_lim'][0], scan_data['scan_info']['scan_lim'][1], num = scan_data['scan_info']['scan_size'])
Show.plot('scanL extended', x_scan, scan_data['scan'])

frequencies = np.arange(4, 8, 1)*1e9
field  = [2500, 2000, 1950, 1900, 1800, 1600, 1400, 1200, 700, 450]
angles  = [  40,   35,   30,   25,   20,   25,   30,   35,  40,  45]
height_field = np.arange(0, len(field))*7.25e7
height_atm   = np.array([   0,    1,  1.2,  1.5,    2,    3,    4,    5,  10,  15], dtype = np.float64)*1e8;
temperature  = np.array([0.01, 0.01,    1,    2,    2,    2,    2,    2,   2,   2], dtype = np.float64)*1e6;
density = 3e15/temperature;

taus = 10**np.linspace(-2, 2, 208)

cost = np.cos(np.deg2rad(angles))
res_los = reo.los_calculate(frequencies
                          , height_field, field, cost
                          , height_atm, temperature, density
                          , taus = taus
                           )

pass


