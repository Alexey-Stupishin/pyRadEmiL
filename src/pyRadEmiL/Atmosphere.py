import ASSolarPy.AtmoProfile as atp
import numpy as np

class Atmosphere:
    NONE = 0
    PLAIN = 1
    MASK = 2
    CUBE = 3

    #-------------------------------------------------------------------------------
    def __init__(self):
        self.__mask = None
        self.__type = self.NONE

    #-------------------------------------------------------------------------------
    @property
    def type(self):
        return self.__type

    #-------------------------------------------------------------------------------
    def atmosphere_mask_set(self, maskN, heights, temperature, density):
        self.__mask[maskN] = atp.AtmoProfile(np.array(heights, dtype = np.float64)
                                           , np.array(temperature, dtype = np.float64)
                                           , np.array(density, dtype = np.float64)
                                            )
        self.__type = self.MASK

    #-------------------------------------------------------------------------------
    def atmosphere_init(self, heights, temperature, density):
        self.__mask = dict()
        self.atmosphere_mask_set('base', heights, temperature, density)
        self.__type = self.PLAIN

    #-------------------------------------------------------------------------------
    def get(self, mask = None):
        if mask is None:
            masksN = None
            Lmask = None
            H, T, D = self.__mask.get('base').get
        else:
            masksN = np.unique(mask)
            nmask = len(masksN)
            N = mask.shape

            maxH = 0
            for value in self.__mask.values():
                maxH = max(maxH, len(value.get[0]))

            size_2D = [nmask, maxH]
            H = np.zeros(size_2D, dtype = np.float64, order="C")
            T = np.zeros(size_2D, dtype = np.float64, order="C")
            D = np.zeros(size_2D, dtype = np.float64, order="C")

            Lmask = np.zeros(nmask, dtype = np.int32, order="C")

            defl = self.__mask.get('base')
            cnt = 0
            for m in masksN:
                Hp, Tp, Dp = self.__mask.get(m, defl).get
                nH = len(Hp)
                Lmask[cnt] = nH
                H[cnt, 0:nH] = Hp
                T[cnt, 0:nH] = Tp
                D[cnt, 0:nH] = Dp
                cnt += 1

            return dict(type = self.__type, mask = mask, masksN = masksN, Lmask = Lmask, H = H, T = T, D = D)
