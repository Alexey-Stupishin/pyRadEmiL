import numpy as np
from numpy import linalg as LA

class Coordinates:
    # [0] - longitude, x
    # [1] - latitude,  y

    #-------------------------------------------------------------------------------
    @staticmethod
    def _lonlat2xy(lonlat):
        y = np.sin(np.deg2rad(lonlat[1]));
        x = np.cos(np.deg2rad(lonlat[1])) * np.sin(np.deg2rad(lonlat[0]));

        return np.array([x, y], dtype = np.float64)

    #-------------------------------------------------------------------------------
    @staticmethod
    def _xy2lonlat(xy):
        latitude = np.rad2deg(np.arcsin(xy[1]))
        longitude = np.rad2deg(np.arcsin(xy[0]/np.sqrt(1-xy[1]**2)))
        # !ToDo check out scope

        return np.array([longitude, latitude], dtype = np.float64)

    #-------------------------------------------------------------------------------
    @staticmethod
    def _dir_cos(lonlat):
        dircos = np.array([0, 0, 0], dtype = np.float64)
        dircos[0] = -np.sin(np.deg2rad(lonlat[0]))
        clon = np.cos(np.deg2rad(lonlat[0]))
        dircos[1] = -clon * np.sin(np.deg2rad(lonlat[1]))
        dircos[2] =  clon * np.cos(np.deg2rad(lonlat[1]))

        return dircos

    #-------------------------------------------------------------------------------
    @staticmethod
    def _rotation_matrix_fov2phs(lonlat):

        sinlon = np.sin(np.deg2rad(lonlat[0]));
        coslon = np.cos(np.deg2rad(lonlat[0]));
        sinlat = np.sin(np.deg2rad(lonlat[1]));
        coslat = np.cos(np.deg2rad(lonlat[1]));

        return np.array([
                         [        coslon,      0,        -sinlon]
                       , [-sinlat*sinlon, coslat, -sinlat*coslon]
                       , [ coslat*sinlon, sinlat,  coslat*coslon]
                        ]
                       , dtype = np.float64 )

    #-------------------------------------------------------------------------------
    def __init__(self, RSun = 960, lonlat = None, xy = None, arc = None, box_coord = None):
        self.__RSun = RSun
        self.__lonlat = lonlat
        self.__xy = xy
        self.__arc = arc

        if not box_coord is None:
            self.__RSun = box_coord['RSUN']
            lat, lon = box_coord['LAT_CEN'], box_coord['LON_CEN']
            self.__lonlat = np.array([lon, lat], dtype = np.float64)
            x, y = box_coord['X_CEN'], box_coord['Y_CEN']
            self.__xy = np.array([x, y], dtype = np.float64)
            # self.__xy = Coordinates._lonlat2xy(self.lonlat)
            self.__arc = self.__xy * self.RSun
        elif not lonlat is None:
            self.__lonlat = np.array(lonlat, dtype = np.float64)
            self.__xy = Coordinates._lonlat2xy(self.__lonlat)
            self.__arc = self.__xy * self.RSun
        elif not xy is None:
            self.__xy = np.array(xy, dtype = np.float64)
            self.__lonlat = Coordinates._xy2lonlat(self.__xy)
            self.__arc = self.__xy * self.RSun
        elif not arc is None:
            self.__arc = np.array(arc, dtype = np.float64)
            self.__xy = self.__arc / self.RSun
            self.__lonlat = Coordinates._xy2lonlat(self.__xy)

    #-------------------------------------------------------------------------------
    @property
    def is_outSun(self):
        # np.sqrt(sum(list(map(lambda _: _**2, self.__xy))))
        return LA.norm(self.__xy) > 1

    #-------------------------------------------------------------------------------
    @property
    def RSun(self):
        return self.__RSun

    #-------------------------------------------------------------------------------
    @property
    def lonlat(self):
        return self.__lonlat
    #-------------------------------------------------------------------------------
    @property
    def xy(self):
        return self.__xy
    #-------------------------------------------------------------------------------
    @property
    def arc(self):
        return self.__arc

    #-------------------------------------------------------------------------------
    @property
    def dir_cos(self):
        return Coordinates._dir_cos(self.__lonlat)

    #-------------------------------------------------------------------------------
    @property
    def rotation_matrix_fov2phs(self):
        return Coordinates._rotation_matrix_fov2phs(self.__lonlat)
