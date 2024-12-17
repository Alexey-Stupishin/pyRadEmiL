import ctypes
import numpy as np
import ASSolarPy.GXBox as gxb
import ASSolarPy.Atmosphere as atm

class REOLibrary:
    #-------------------------------------------------------------------------------
    @staticmethod
    def _intensity2temp(intensity, frequency):
        return 6.513e36 * intensity/frequency**2

    #-------------------------------------------------------------------------------
    @staticmethod
    def _fluxpixel2temp(flux, frequency, step = 1):
        return REOLibrary._intensity2temp(flux/(2.35e8*step**2), frequency)

    #-------------------------------------------------------------------------------
    @staticmethod
    def _temp2intensity(temperature, frequency):
        return 1.535e-37 * temperature*frequency**2

    #-------------------------------------------------------------------------------
    @staticmethod
    def _temp2fluxpixel(temperature, frequency, step = 1):
        return REOLibrary._temp2intensity(temperature, frequency) * 2.35e8 * step**2

    #-------------------------------------------------------------------------------
    def __init__(self, lib_path):
        self.__mptr1 = np.ctypeslib.ndpointer(dtype = np.float64, ndim = 1, flags = "C")
        self.__mptr2 = np.ctypeslib.ndpointer(dtype = np.float64, ndim = 2, flags = "C")
        self.__mptr3 = np.ctypeslib.ndpointer(dtype = np.float64, ndim = 3, flags = "C")
        self.__mpint1 = np.ctypeslib.ndpointer(dtype = np.int32, ndim = 1, flags = "C")
        self.__mpint2 = np.ctypeslib.ndpointer(dtype = np.int32, ndim = 2, flags = "C")
        self.__mpint3 = np.ctypeslib.ndpointer(dtype = np.int32, ndim = 3, flags = "C")
        self.__mpstr = ctypes.POINTER(ctypes.c_char)
        self.__mvoid = ctypes.c_void_p
        self.__mint = ctypes.c_int32
        self.__mreal = ctypes.c_double
        self.__mdw = ctypes.c_uint32
        self.__mpdw1 = np.ctypeslib.ndpointer(dtype = np.uint32, ndim = 1, flags = "C")

        self.__map_size = np.array([0, 0], dtype = np.int32, order="C")
 
        libc_reo = ctypes.CDLL(lib_path)

        create_func = libc_reo.utilInitialize
        create_func.argtypes = [self.__mvoid, self.__mvoid]
        create_func.restype = self.__mdw

        set_int_func = libc_reo.utilSetInt
        set_int_func.argtypes = [self.__mpstr, self.__mint]
        set_int_func.restype = self.__mint

        set_double_func = libc_reo.utilSetDouble
        set_double_func.argtypes = [self.__mpstr, self.__mreal]
        set_double_func.restype = self.__mint

        field_set_func = libc_reo.physSetField
        field_set_func.argtypes = [self.__mptr3, self.__mptr3, self.__mptr3, self.__mvoid, self.__mvoid, self.__mvoid, self.__mvoid, self.__mvoid]
        field_set_func.restype = self.__mdw

        field_markup_func = libc_reo.physSetFieldPlaneMeshgridRect
        field_markup_func.argtypes = [self.__mvoid, self.__mreal, self.__mpint1, self.__mptr1, self.__mint]
        field_markup_func.restype = self.__mdw
     
        atmosphere_plain_func = libc_reo.physSetAtmosphere
        atmosphere_plain_func.argtypes = [self.__mint, self.__mptr1, self.__mptr1, self.__mptr1]
        atmosphere_plain_func.restype = self.__mdw
     
        atmosphere_mask_func = libc_reo.physSetAtmosphereMask
        atmosphere_mask_func.argtypes = [self.__mint, self.__mint, self.__mpint1, self.__mptr2, self.__mptr2, self.__mptr2, self.__mpint1, self.__mpint1, self.__mpint2]
        atmosphere_mask_func.restype = self.__mdw
     
        map_calculate_func = libc_reo.gstcCalculateMapAll
        # map_calculate_func.argtypes sets dynamically
        map_calculate_func.restype = self.__mdw

        los_calculate_func = libc_reo.gstcCalculateLineFromLOS
        los_calculate_func.argtypes = [self.__mint, self.__mptr1, self.__mptr1, self.__mptr1 # int L, double *H /* L */, double *B, double *cost,
                                       , self.__mint, self.__mptr1, self.__mptr1, self.__mptr1 # int La, double *Ha, double *_T, double *_N,
                                       , self.__mint, self.__mptr1, self.__mint, self.__mpint1 # int nf, double *pf, int ns4calc, int *s4calc,
                                       , self.__mint, self.__mptr1                             # int nLayers, double *pTauL /* nLayers */,
                                       , self.__mpint2, self.__mptr2, self.__mptr2             #  int *depth /* 2*nf */, double *pF /* full intensity, 2*nf */, double *pTau /* full tau, 2*nf */,
                                       , self.__mptr3, self.__mptr3, self.__mpint3             #  double *pHL /* heights, 2*nLayers*nf */, double *pFL /* intensity, 2*nLayers*nf */, int *psL /* 2*nLayers*nf */
                                       , self.__mpdw1                                          # DWORD *pRC
                                       ]
        los_calculate_func.restype = self.__mdw

        map_get_scan_limits_func = libc_reo.gstcGetScanLimits
        map_get_scan_limits_func.argtypes = [self.__mptr1, self.__mpint1]
        map_get_scan_limits_func.restype = self.__mdw

        map_convolve_func = libc_reo.gstcConvolve
        # map_convolve_func.argtypes sets dynamically
        map_convolve_func.restype = self.__mdw

        get_version_func = libc_reo.utilGetVersion
        get_version_func.argtypes = [self.__mpstr, self.__mint]
        get_version_func.restype = self.__mint

        self.__lib_path = lib_path
        self.__func_set = {'create_func':create_func
                         , 'set_int_func':set_int_func
                         , 'set_double_func':set_double_func
                         , 'field_set_func':field_set_func
                         , 'field_markup_func':field_markup_func
                         , 'atmosphere_plain_func':atmosphere_plain_func
                         , 'atmosphere_mask_func':atmosphere_mask_func
                         , 'map_calculate_func':map_calculate_func
                         , 'los_calculate_func':los_calculate_func
                         , 'map_get_scan_limits_func':map_get_scan_limits_func
                         , 'map_convolve_func':map_convolve_func
                         , 'get_version_func':get_version_func
                          }
        self.__pointer = create_func(0, 0)
        self.__atmo_mask = dict()
        self.__arc_box = None

    #-------------------------------------------------------------------------------
    def set_int(self, prop, vint):
        return self.__func_set['set_int_func'](prop.encode('utf-8'), vint)

    #-------------------------------------------------------------------------------
    def set_double(self, prop, vdouble):
        return self.__func_set['set_double_func'](prop.encode('utf-8'), vdouble)

    #-------------------------------------------------------------------------------
    @property
    def get_version(self):
        buflen = 512
        bufenc = ''.ljust(buflen).encode('utf-8')
        rc = self.__func_set['get_version_func'](bufenc, buflen)
        buffer = bufenc.decode('utf-8')
        term = buffer.find(chr(0))

        return buffer[0:term]

    #-------------------------------------------------------------------------------
    def __phys_field_set(self, box):
        # raise box is of gxb class
        cube = box.get_cube(order = 'C')
        Nc = cube['bx'].shape
        N = np.array([Nc[2], Nc[1], Nc[0]], dtype = np.int32)
        vcosc = box.get_index_value('VCOS')
        vcos = np.array([vcosc[1], vcosc[0], vcosc[2]], dtype = np.float64)
        dx = box.get_index_value('DY')
        dy = box.get_index_value('DX')
        mod_step = np.array([dx, dy, dx], dtype = np.float64)
        mod_base = - (N[0:2]-1)/2*mod_step[0:2]
 
        return self.__func_set['field_set_func'](cube['bx'], cube['by'], cube['bz'], N.tobytes(), N.tobytes(), 
                                                 vcos.tobytes(), mod_step.tobytes(), mod_base.tobytes())

    #-------------------------------------------------------------------------------
    def __phys_field_markup(self, R_arc, map_step_arc, pos_angle0, setFOV = 0):
        map_stepR = map_step_arc/R_arc
        map_step = np.array([map_stepR, map_stepR], dtype = np.float64)
        pos_angle = float(pos_angle0)
        map_size = np.array([0, 0], dtype = np.int32, order="C")
        map_base = np.array([0, 0], dtype = np.float64, order="C")

        rc = self.__func_set['field_markup_func'](map_step.tobytes(), pos_angle, map_size, map_base, int(not setFOV))

        self.__map_step = map_step_arc
        self.__map_size = map_size
        self.__map_base = map_base
        arc_box = np.array([[map_base[0]*R_arc, map_base[0]*R_arc+(map_size[0]-1)*map_step_arc]
                          , [map_base[1]*R_arc, map_base[1]*R_arc+(map_size[1]-1)*map_step_arc]]
                          , dtype = np.float64, order="C")
        self.__arc_box = arc_box
    
        # field = 'reoGetField'
        return dict(map_size = map_size[::-1]
                  , map_base = map_base[::-1]
                  , R_arc = R_arc
                  , arc_box = arc_box[::-1, :]
                  , rc = rc)

    #-------------------------------------------------------------------------------
    def prepare(self, box, map_step, pos_angle
                           #  , setFOV = 0
                           #  , version_info = version_info
                           , freefree = 1
                           , distribution_type = 1
                           , distribution_kappaK = 5.0
                           , use_LaplasMethod = 1
                           , use_QT = 1
                             ):

        self.set_int('cycloCalc.ConsiderFreeFree', freefree)
        self.set_int('cycloCalc.Distribution.Type', distribution_type)
        self.set_double('cycloCalc.Distribution.kappaK', distribution_kappaK)
        self.set_int('cycloCalc.LaplasMethod.Use', use_LaplasMethod)
        self.set_int('cycloLine.ZoneSearch.QT.Use', use_QT)

        setFOV = 0

        self.__phys_field_set(box)
        R_arc = box.get_index_value('RSUN')
        return self.__phys_field_markup(R_arc, map_step, pos_angle, setFOV)

    #-------------------------------------------------------------------------------
    def __def_scan_limits(self, scan_lim):
        if scan_lim is None:
            scan_lim = self.__arc_box[1, :]
            scan_pos = [0, self.__map_size[1]-1]
            scan_pos_pass = 0
            type_scan_pos = self.__mvoid
            scan_size = self.__map_size[1]
        else:
            #assert self.__arc_box is None

            scan_lim = np.array(scan_lim, dtype = np.float64)
            scan_pos = np.array([0, 0], dtype = np.int32)
            type_scan_pos = self.__mpint1
            rc = self.__func_set['map_get_scan_limits_func'](scan_lim, scan_pos)
            scan_pos_pass = scan_pos
            scan_size = scan_pos[1] - scan_pos[0] + 1

        return dict(scan_lim = scan_lim
                  , scan_pos = scan_pos
                  , scan_size = scan_size
                  , scan_pos_pass = scan_pos_pass
                  , type_scan_pos = type_scan_pos
                   )

    #-------------------------------------------------------------------------------
    def map_calculate(self, frequency, harmonics = [2, 3, 4], taus = [100], view_mask = None
                         , mode = 6, n_beam = 1, beam_c = 0.0, beam_b = 0.0, scan_lim = None, map_unit = 'sfu'):

        # assert self.__arc_box is None

        harmonics = np.array(harmonics, dtype = np.int32)
        nharm = len(harmonics) 
        taus = np.array(taus, dtype = np.float64)
        ntaus = len(taus)

        map_calculate_func = self.__func_set['map_calculate_func']
        type_view_mask = self.__mpint2
        if view_mask is None:
            view_mask = 0
            type_view_mask = self.__mvoid

        sp = self.__def_scan_limits(scan_lim)

        map_calculate_func.argtypes = [self.__mreal, self.__mint, self.__mpint1, self.__mint, self.__mptr1, type_view_mask
                                     , self.__mint, self.__mint, self.__mreal, self.__mreal
                                     , self.__mint, self.__mint
                                     , self.__mpint2, self.__mptr2, self.__mptr2, self.__mptr3, self.__mptr3, self.__mpint3
                                     , self.__mpint2, self.__mptr2, self.__mptr2, self.__mptr3, self.__mptr3, self.__mpint3
                                     , self.__mptr1, self.__mptr1, sp['type_scan_pos']
                                     , self.__mvoid, self.__mvoid, self.__mvoid, self.__mvoid, self.__mvoid]
        
        size_2D = np.flip(self.__map_size)
        size_3D = np.flip(np.append(self.__map_size, ntaus))
        depthR = np.zeros(size_2D, dtype = np.int32, order="C")
        fluxR = np.zeros(size_2D, dtype = np.float64, order="C")
        tauR = np.zeros(size_2D, dtype = np.float64, order="C")
        heightsR = np.zeros(size_3D, dtype = np.float64, order="C")
        h_fluxR = np.zeros(size_3D, dtype = np.float64, order="C")
        h_harmR = np.zeros(size_3D, dtype = np.int32, order="C")
        depthL = np.zeros(size_2D, dtype = np.int32, order="C")
        fluxL = np.zeros(size_2D, dtype = np.float64, order="C")
        tauL = np.zeros(size_2D, dtype = np.float64, order="C")
        heightsL = np.zeros(size_3D, dtype = np.float64, order="C")
        h_fluxL = np.zeros(size_3D, dtype = np.float64, order="C")
        h_harmL = np.zeros(size_3D, dtype = np.int32, order="C")
        scanR = np.zeros(sp['scan_size'], dtype = np.float64, order="C")
        scanL = np.zeros(sp['scan_size'], dtype = np.float64, order="C")

        rc = map_calculate_func(frequency, nharm, harmonics, ntaus, taus, view_mask, mode, n_beam, beam_c, beam_b
                                                 , -1, -1
                                                 , depthR, fluxR, tauR, heightsR, h_fluxR, h_harmR
                                                 , depthL, fluxL, tauL, heightsL, h_fluxL, h_harmL
                                                 , scanR, scanL, sp['scan_pos_pass']
                                                 , 0, 0, 0, 0, 0
                                                  )

        if map_unit == 'T':
            fluxR = self.fluxpixel2temp(fluxR, frequency)
            fluxL = self.fluxpixel2temp(fluxL, frequency)
            h_fluxR = self.fluxpixel2temp(h_fluxR, frequency)
            h_fluxL = self.fluxpixel2temp(h_fluxL, frequency)

        return dict(depthR = depthR.transpose()
                  , fluxR = fluxR.transpose()
                  , tauR = tauR.transpose() 
                  , heightsR = heightsR.transpose()
                  , h_fluxR = h_fluxR.transpose()
                  , h_harmR = h_harmR.transpose()
                  , depthL = depthL.transpose()
                  , fluxL = fluxL.transpose()
                  , tauL = tauL.transpose()
                  , heightsL = heightsL.transpose()
                  , h_fluxL = h_fluxL.transpose()
                  , h_harmL = h_harmL.transpose()
                  , scanR = scanR
                  , scanL = scanL
                  , scan_info = sp
                  , rc = rc
                   )
    
    #-------------------------------------------------------------------------------
    def set_atmosphere(self, atm, mask = None):
        # raise: atmo should be an Atmosphere (atm)
        atmo = atm.get(mask)
        # raise: Atmosphere should have at least 'base'
        H = atmo['H']
        T = atmo['T']
        D = atmo['D']
        if atm.type == atm.PLAIN:
            return self.__func_set['atmosphere_plain_func'](len(H), H, T, D)
        else:
            masksN = atmo['masksN']
            Lmask = atmo['Lmask']

            mask = atmo['mask']
            mask = np.transpose(mask).astype(np.int32, order="C")

            par_size = H.shape
            mask_size = mask.shape

            return self.__func_set['atmosphere_mask_func'](par_size[0], par_size[1], Lmask, H, T, D, masksN, np.asarray(mask_size, dtype = np.int32, order="C"), mask)

    #-----------------------------------------------------------------------------------------
    def map_convolve(self, rmap, frequency, mode = 6, n_beam = 1, beam_c = 0.0, beam_b = 0.0, scan_lim = None):
        sp = self.__def_scan_limits(scan_lim)

        self.__func_set['map_convolve_func'].argtypes = [self.__mptr2, self.__mreal, self.__mptr1
                                                       , self.__mint, self.__mint, self.__mreal, self.__mreal
                                                       , sp['type_scan_pos']
                                                        ]

        scan = np.zeros(sp['scan_size'], dtype = np.float64, order="C")
        rc = self.__func_set['map_convolve_func'](rmap.transpose().astype(np.float64, order="C"), frequency, scan, mode, n_beam, beam_c, beam_b, sp['scan_pos'])

        return dict(scan = scan, scan_info = sp)

    #-------------------------------------------------------------------------------
    def fluxpixel2temp(self, flux, frequency):
        return self._fluxpixel2temp(flux, frequency, self.__map_step)

    #-------------------------------------------------------------------------------
    def los_calculate(self, frequencies
                    , heights_field, field, cost
                    , heights_atm, temperature, density
                    , harmonics = [2, 3, 4], taus = [100]):
        n_freq = len(frequencies)
        n_taus = len(taus)
        #size_2D = np.array([2, n_freq], dtype = np.int32, order="C")
        #size_3D = np.array([2, n_taus, n_freq], dtype = np.int32, order="C")
        size_2D = np.array([n_freq, 2], dtype = np.int32, order="C")
        size_3D = np.array([n_freq, n_taus, 2], dtype = np.int32, order="C")
        depth = np.zeros(size_2D, dtype = np.int32, order="C")
        flux = np.zeros(size_2D, dtype = np.float64, order="C")
        tau = np.zeros(size_2D, dtype = np.float64, order="C")
        heights = np.zeros(size_3D, dtype = np.float64, order="C")
        h_flux = np.zeros(size_3D, dtype = np.float64, order="C")
        h_harm = np.zeros(size_3D, dtype = np.int32, order="C")

        n_harm = len(harmonics)
        n_field = len(heights_field) # should be == len(field) == len(cost), raise it
        n_atm = len(heights_atm) # should be == len(temperature) == len(density), raise it
        rc = np.array([1], dtype = np.uint32);
        rc = self.__func_set['los_calculate_func'](
                  n_field, np.array(heights_field, dtype = np.float64), np.array(field, dtype = np.float64), np.array(cost, dtype = np.float64)
                , n_atm, np.array(heights_atm, dtype = np.float64), np.array(temperature, dtype = np.float64), np.array(density, dtype = np.float64)
                , n_freq, np.array(frequencies, dtype = np.float64), n_harm, np.array(harmonics, dtype = np.int32)
                , n_taus, np.array(taus, dtype = np.float64)
                , depth, flux, tau
                , heights, h_flux, h_harm
                , rc
                                                  )            
        return dict(depth = depth, flux = flux, tau = tau, heights = heights, h_flux = h_flux, h_harm = h_harm)

