import numpy as np
import numpy as np

import ASSolarPy.AtmoProfile as atp
import ASSolarPy.REOwrapper as REOwrapper
import ASSolarPy.AtmoProfile as atp

class QuietSun:
    Zirin1991 = 1

    _Zirin1991 = np.array([ 
                [1.4, 	 70.5]
              , [1.6, 	 63.8]
              , [1.8, 	 52.2]
              , [2.0, 	 42.9]
              , [2.4, 	 32.8]
              , [2.8, 	 27.1]
              , [3.2, 	 24.2]
              , [3.6, 	 21.7] 
              , [4.2, 	 19.4] 
              , [5.0, 	 17.6] 
              , [5.8, 	 15.9] 
              , [7.0, 	 14.1] 
              , [8.2, 	 12.9] 
              , [9.4, 	 12.2] 
              , [10.6, 	 11.3] 
              , [11.8,	 11.0] 
              , [13.2,	 10.8] 
              , [14.8,	 10.8] 
              , [16.4,	 10.7]
              , [18.0, 	 10.3]
            ])

    #-------------------------------------------------------------------------------
    def __init__(self):
        pass

    #-------------------------------------------------------------------------------
    @staticmethod
    def _QS_temperature(frequency, model = Zirin1991):
        return np.interp(frequency*1e-9, QuietSun._Zirin1991[:, 0], QuietSun._Zirin1991[:, 1]) * 1000

    #-------------------------------------------------------------------------------
    @staticmethod
    def _map_flat(pos, intensity = 1):
        xg, yg = np.meshgrid(pos, pos)
        dist = np.sqrt(xg**2 + yg**2)
        flux_map = np.ones_like(dist)*intensity
        flux_map[dist > 1] = 0

        return flux_map

        pass

    #-------------------------------------------------------------------------------
    def calc_freefree_cut(self, profile, reo, posx, frequencies):
        height_field = [1]
        field = [0]
        cost = [0]

        reo.set_int('cycloCalc.ConsiderFreeFree', 1)
        reo.set_int('cycloCalc.ConsiderFreeFree.Only', 1)

        flux = np.zeros([len(frequencies), len(posx)])
        for kx in range(len(posx)):
            height_atm, temperature, density = profile.los(posx[kx], 0)
            if not height_atm is None:
                res_los = reo.los_calculate(frequencies
                              , height_field, field, cost
                              , height_atm, temperature, density
                               )

                flux[:, kx] = res_los['flux'][:, 0]

        return flux

    #-------------------------------------------------------------------------------
    def calc_freefree_map(self, profile, reo, step, max_R, frequencies):
        posx = np.arange(0, max_R+step, step)
        flux00 = self.calc_freefree_cut(profile, reo, posx, frequencies)
        flux0 = np.hstack([np.flip(flux00, 1), flux00[:, 1::]])
        map_size = 2*len(posx) - 1
        map_center = len(posx) - 1

        flux_map = np.zeros([len(frequencies), map_size, map_size])
        flux_map[:, map_center, :] = flux0

        pos = np.hstack([- np.flip(posx), posx[1::]])
        xg, yg = np.meshgrid(pos, pos)
        dist = np.sqrt(xg**2 + yg**2)

        for kf in range(len(frequencies)):
            flux_map[kf,:,:] = np.interp(dist, pos, flux0[kf, :], right = 0)

        return (flux_map, pos)
        