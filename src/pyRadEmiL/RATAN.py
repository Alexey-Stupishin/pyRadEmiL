from astropy.io import fits
import numpy as np
import scipy

class RATAN:
    FILE_FITS = 0
    FILE_SKAT = 1

    BEAM_SIMPLE  = 1 # simple, from SAO site
    BEAM_NARROW  = 2 # theoretical for spiral feed, unrealistic
    BEAM_STOH    = 3 # some last private letter (S. Tokchukova? Need to clarify)
    BEAM_UDEF    = 4 # user defined by table lookup, TODO
    BEAM_WIDE    = 5 # wide, source unknown
    BEAM_2023    = 6 # 1-3 GHz, 2023 (spring) by N.Ovchinnikova/M.Lebedev

    #-------------------------------------------------------------------------------
    @staticmethod
    def _position_angle(azimuth, sol_dec, solar_p):
        if solar_p > 180:
           solar_p -= 360
        if solar_p < -180:
           solar_p += 360

        return solar_p + np.rad2deg(np.arcsin(-np.tan(np.deg2rad(azimuth))* np.tan(np.deg2rad(sol_dec))))

    #-------------------------------------------------------------------------------
    @staticmethod
    def beam(frequencies, mode = BEAM_2023, c = 0, b = 0):
        wave = 30*1e9/frequencies
        vert = 7.5*60*wave

        match mode:
            case RATAN.BEAM_SIMPLE:
                horz =          8.5  *wave
            case RATAN.BEAM_NARROW:
                horz =  4.38  + 6.87 *wave
            case RATAN.BEAM_STOH:
                horz =  0.009 + 8.338*wave
            case RATAN.BEAM_UDEF:
                horz = 0 # TODO
            case RATAN.BEAM_WIDE:
                horz =  0.2   + 9.4  *wave
            case RATAN.BEAM_2023:
                horz = -0.16  + 8.162*wave
            case _:
                horz = 0 # raise?

        return dict(horz = horz, vert = vert)

    #-------------------------------------------------------------------------------
    @staticmethod
    def __gauss_norm(halfwidth, pos, length):
        m = np.log(2)/halfwidth**2
        k = np.arange(0, length)

        return np.sqrt(m/np.pi)*np.exp(-((k-pos)**2)*m)

    #-------------------------------------------------------------------------------
    @staticmethod
    def gauss_test(halfwidth, pos, length):
        sigma = halfwidth*2;
        x = np.arange(0, length) - pos;
        gauss = np.exp(-((x / (sigma / (2 * np.sqrt(2 * np.log(2))))) ** 2) / 2) / (
            sigma / (2 * np.sqrt(2 * np.log(2)))  * np.sqrt(2 * np.pi) 
                     )
        return gauss/np.sum(gauss)
    
    #-------------------------------------------------------------------------------
    @staticmethod
    def twosinc1(freqGHz, pos, length):
        x = np.arange(0, length) - pos;
        xx = x / 2 / 3600 / 180 * np.pi
        k = 2 * np.pi * freqGHz / scipy.constants.speed_of_light
        y = k * xx
        a = 300
        b = a / (2 * np.arccos(.2))
        return freqGHz * (2*b*(-np.cos((a*y)/2)*np.sin(a/(2*b)) + b*y*np.cos(a/(2*b))*np.sin((a*y)/2)))/(-1 + b**2*y**2)

    #-------------------------------------------------------------------------------
    @staticmethod
    def twosinc_sq1(freqGHz, pos, length):
        sinc = np.abs(RATAN.twosinc1(freqGHz, pos, length)) ** 2
        return sinc / np.sum(sinc)

    #-------------------------------------------------------------------------------
    @staticmethod
    def __gauss_points(halfwidth, points):

        a0 = points[0];
        n = len(points) - 1;
        P = (points[n] - a0)/n;
        m = np.log(2)/halfwidth**2
        k = np.arange(0, len(points))

        return np.exp(-((k*P + a0)**2)*m)

    #-------------------------------------------------------------------------------
    @staticmethod
    def __beams(frequency, sizes, step, base, mode = BEAM_2023, c = 0, b = 0):
        length = sizes[0]
        points = base[1]/step[1] + np.arange(0, sizes[1])
        hpbw = RATAN.beam(frequency, mode, c, b)
        beamH = RATAN.__gauss_norm(hpbw['horz']/step[0]/2, (length-1)/2 + length, 3*length)
        beamH1 = RATAN.gauss_test(hpbw['horz']/step[0]/2, (length-1)/2 + length, 3*length)
        beamH2 = RATAN.twosinc_sq1(frequency, (length-1)/2 + length, 3*length)
        ss = 4095
        beamH2x = RATAN.twosinc_sq1(frequency, ss//2, ss)
        beamV = RATAN.__gauss_points(hpbw['vert']/step[1]/2, points)

        return dict(beamH = beamH, beamV = beamV)
    
    #-------------------------------------------------------------------------------
    def __init__(self, file, type = FILE_FITS):
        hdul = fits.open(file)
        self.__header = hdul[0].header
        self.__data = hdul[0].data
        pass

    #-------------------------------------------------------------------------------
    def get_index_value(self, key):
        return self.__header[key]

    #-------------------------------------------------------------------------------
    @property
    def position_angle(self):
        return self._position_angle(self.get_index_value('AZIMUTH'), self.get_index_value('SOL_DEC'), self.get_index_value('SOLAR_P'))

    #-------------------------------------------------------------------------------
    @staticmethod
    def _convolve(flux_map, frequency, vis_step, arc_base, mode = BEAM_2023, c = 0, b = 0):
        sizes = np.flip(flux_map.shape)
        beams =  RATAN.__beams(frequency, sizes, vis_step, [0, arc_base[0]], mode, c, b) # arc_base[1] ?
        beamH = beams['beamH']
        beamV = beams['beamV']

        pure_scan = flux_map.T.dot(beamV)

        n_horz = len(beamH)
        n_scan = len(pure_scan)
        nzeros = np.floor((n_horz-n_scan)/2).astype(np.int32)
        nlast = n_horz - n_scan

        pure_scan_ex = np.hstack([np.zeros(nzeros), pure_scan, np.zeros(n_horz - nzeros - n_scan)])

        scan_ex = np.convolve(pure_scan_ex, beamH, mode = 'same')
        scan = scan_ex[nzeros:nlast] / vis_step[0] # vis_step[1] ?

        return scan

        pass
