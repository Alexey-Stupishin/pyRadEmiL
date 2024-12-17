from scipy.io import readsav
import numpy as np
import ASSolarPy.Common as asc
from ASSolarPy.Coordinates import Coordinates 

class GXBox(Coordinates):
    #-------------------------------------------------------------------------------
    @staticmethod
    def _as_dict(record_):
        return {name:record_[name] for name in record_.dtype.names}

    #-------------------------------------------------------------------------------
    @staticmethod
    def _get_details_mask(bz, cont):
        abz = np.abs(bz)

        mag_qs = 10.0  # 10 Gauss for QS
        thr_plage = 3 # MF in plage is thr_plage times stronger than QS

        qs = abz < mag_qs
        cutoff_qs = np.sum(cont[qs])/np.count_nonzero(qs)

        # exclude sunspots
        non_zero = cont > 0
        n_used = np.count_nonzero(non_zero)
        sub = np.logical_and(cont > cutoff_qs*0.9, non_zero)
        sub_count = np.count_nonzero(sub)
        hist, edges = np.histogram(cont[sub], bins = sub_count)
        hist_sum = np.cumsum(hist)/sub_count
        idx_b = np.argmin(np.abs(hist_sum - 0.75))
        cutoff_b = edges[idx_b];
        idx_f = np.argmin(np.abs(hist_sum - 0.97))
        cutoff_f = edges[idx_f];

        model_mask = np.zeros(cont.shape, dtype = np.int32, order="C")

        model_mask[cont <= 0.65*cutoff_qs] = 7 # umbra
        model_mask[np.logical_and(cont > 0.65*cutoff_qs, cont <= 0.9*cutoff_qs)] = 6 # penumbra
        model_mask[np.logical_and(cont > cutoff_f, cont <= 1.19*cutoff_qs)] = 3 # enhanced NW
        model_mask[np.logical_and(cont > cutoff_b, cont <= cutoff_f)] = 2 # NW lane
        model_mask[np.logical_and(cont > 0.9*cutoff_qs, cont <= cutoff_b)] = 1 # IN
        model_mask[np.logical_and(cont > 0.95*cutoff_qs, np.logical_and(cont <= cutoff_f, abz > thr_plage*mag_qs))] = 4 # plage
        model_mask[np.logical_and(cont > 1.01*cutoff_qs, abz > thr_plage*mag_qs)] = 5 # facula

        return model_mask    

    #-------------------------------------------------------------------------------
    def __init__(self, filename):
        sav_data = readsav(filename)
        # keys = sav_data.keys()
        box = sav_data['box']
        self.__box = self._as_dict(box[0])

        RSUN = self.get_index_value('RSUN')
        LAT_CEN, LON_CEN = self.get_index_value('LAT_CEN'), self.get_index_value('LON_CEN')
        X_CEN, Y_CEN = self.get_index_value('X_CEN'), self.get_index_value('Y_CEN')
        box_coord = dict(RSUN = RSUN, LAT_CEN = LAT_CEN, LON_CEN = LON_CEN, X_CEN = X_CEN, Y_CEN = Y_CEN)
        super().__init__(box_coord = box_coord)

    #-------------------------------------------------------------------------------
    def get_cube(self, order = 'F'):
        cube = {'bx':self.__box['BX'], 'by':self.__box['BY'], 'bz':self.__box['BZ']}
        if order == 'C':
            by = np.transpose(self.__box['BX'], (0, 2, 1)).astype(np.float64, order="C")
            bx = np.transpose(self.__box['BY'], (0, 2, 1)).astype(np.float64, order="C")
            bz = np.transpose(self.__box['BZ'], (0, 2, 1)).astype(np.float64, order="C")
            cube = {'bx':bx, 'by':by, 'bz':bz}

        return cube

    #-------------------------------------------------------------------------------
    def get_index_value(self, key):
        index = self.__box['INDEX'][0]
        return index[key]

    #-------------------------------------------------------------------------------
    @property
    def get_details_mask(self):
        base = self._as_dict(self.__box['BASE'][0])
        bz = base['BZ'].astype(np.float64)
        cont = base['IC'].astype(np.float64)

        return self._get_details_mask(bz, cont)
