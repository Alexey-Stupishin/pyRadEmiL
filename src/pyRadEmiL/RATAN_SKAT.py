import re
import numpy as np

class RATAN_SKAT:
    MODE_NONE = 0
    MODE_SCANS = 1
    MODE_SPECTRA = 2
    
    #-------------------------------------------------------------------------------
    @staticmethod
    def __convert_head(key, value):
        as_string = ('version', 'created', 'user_sel_freqs', 'DATE-OBS', 'TIME-OBS', 'CAL_BY', 'METHOD', 'PTS_BY', 'UNITS', '')        
        as_int = ('N_FREQS', 'N_POINTS', 'N_POS')        
        if key not in as_string:
            if key in as_int:
                value = int(value)
            else:
                value = float(value)
            
        return value
            
    #-------------------------------------------------------------------------------
    def __read_head(self, line):
        sect = 'none'
        key = ''
        value = 0
        mode = self.MODE_NONE
    
        predefined = {'RATAN Spectra-At-Positions Data File v ':'version',
                      'RATAN Selected Scans Data File v ':'version',
                      'generated at ':'created',
                      'frequencies selected by user for scan':'user_sel_freqs',
                     }
        
        keys = list(predefined.keys())
        for predef in keys:
            pattern = '# +' + predef + '(.*)'
            found = re.match(pattern, line)
            if found is not None:
                sect = 'nonkey'
                key = predefined[predef]
                value = found.group(1)
                mode = self.MODE_SPECTRA if keys.index(predef) == 0 else \
                      (self.MODE_SCANS if keys.index(predef) == 1 else self.MODE_NONE)
                return {'sect':sect, 'key':key, 'value':RATAN_SKAT.__convert_head(key, value), 'mode':mode}
 
        pattern = '# +PAR.(\S*) *= *(\S*)'
        found = re.match(pattern, line)
        if found is not None:
            sect = 'par'
            key = found.group(1)
            value = found.group(2)
            return {'sect':sect, 'key':key, 'value':RATAN_SKAT.__convert_head(key, value)}
            
        pattern = '# +(\S*) *= *(\S*)'
        found = re.match(pattern, line)
        if found is not None:
            sect = 'main'
            key = found.group(1)
            value = found.group(2)
            return {'sect':sect, 'key':key, 'value':RATAN_SKAT.__convert_head(key, value)}

        return {'sect':sect, 'key':key, 'value':RATAN_SKAT.__convert_head(key, value), 'mode':mode}
        
    #-------------------------------------------------------------------------------
    def __init__(self):
        self._mode = self.MODE_NONE
        self._header = {}
        self._params = {}
        self._n_pos = None
        self._pos = None
        self._n_freqs = None
        self._freqs = None
        self._right = None
        self._left = None
        self._qs = None
        self._ffpos = None
        self._ff = None

    #-------------------------------------------------------------------------------
    def load(self, filename):
        f = open(filename, 'r') # ToDo: if not exist?
    
        cnt = 0        

        fstate = 'header'
        skip = False
        while True:
            if not skip:
                fstr = f.readline()
                if len(fstr) == 0:
                    break
                
            skip = False
            match fstate:
                case 'header':
                    result = self.__read_head(fstr)
                    match result['sect']:
                        case 'none':
                            fstate = 'headinfo'
                            skip = True
                        case 'nonkey':
                            if self._mode == self.MODE_NONE:
                                self._mode = result['mode']
                            self._header[result['key']] = result['value']
                        case 'main':
                            key = result['key']
                            if key == 'N_POINTS':
                                key = 'N_POS'
                            self._header[key] = result['value']
                        case 'par':
                            self._params[result['key']] = result['value']
                
                case _:
                    s = list(map(float, fstr.split()))
                    if fstate == 'headinfo':
                        self._n_pos = self._header['N_POS']
                        self._n_freqs = self._header['N_FREQS']

                        self._right = np.zeros((self._n_pos, self._n_freqs))
                        self._left  = np.zeros((self._n_pos, self._n_freqs))
                        self._qs    = np.zeros((self._n_pos, self._n_freqs))
                            
                        par = np.array(s)
                        
                        if self._mode == self.MODE_SCANS:
                            self._freqs = par[1::3]
                            self._pos = np.zeros((1, self._n_pos))
                        else:
                            self._freqs = np.zeros(self._n_freqs)
                            self._ffpos = par[-1]
                            par = par[:-2]
                            self._pos = par[1::3]
                            self._ff = np.zeros((1, self._n_freqs))

                        fstate = 'table'
                    else: # table
                        vstr = np.array(s)
                        if self._mode == self.MODE_SCANS:
                            self._pos[cnt] = vstr[0]
                            self._right[cnt, :] = vstr[1::3]
                            self._left[cnt, :] = vstr[2::3]
                            self._qs[cnt, :] = vstr[3::3]
                        else:
                            self._freqs[cnt] = vstr[0]
                            self._ff = vstr[-1]
                            vstr = vstr[:-1]
                            self._right[:, cnt] = np.reshape(vstr[1::3], (1, self._n_pos))
                            self._left[:, cnt] = np.reshape(vstr[2::3], (1, self._n_pos))
                            self._qs[:, cnt] = np.reshape(vstr[3::3], (1, self._n_pos))
                        
                        cnt += 1

        pass

    # ratan = 
    #   mode: 'spectra'
    # header: [1x1 struct]
    # params: [1x1 struct]
    #    pos: [9x1 double]
    #  freqs: [1x71 double]
    #  right: [9x71 double]
    #   left: [9x71 double]
    #     qs: [9x71 double]
    #  ffpos: 131.09534
    #     ff: [1x71 double]
    #  version: 'v 1.0.20.52'
    #  created: 'Wed Sep 02 11:48:40 202'
    # DATE_OBS: '2015/12/18'
    # TIME_OBS: '08:26:24.100'
    #   CDELT1: 2.8812163
    #  AZIMUTH: 10
    #  SOL_DEC: -23.379
    #  SOLAR_R: 975.25
    #  SOLAR_P: -351.10001
    #  SOLAR_B: -1.3
    #  RATAN_P: 13.271707
    #    N_POS: 9
    #  N_FREQS: 71
    #    shift: 0

    #   mode: 'scans'
    # header: [1x1 struct]
    # params: [1x1 struct]
    #    pos: [234x1 double]
    #  freqs: [3000000000 4000000000 6000000000 7900000100 10000000000 11900000000 14000000000 16500000000]
    #  right: [234x8 double]
    #   left: [234x8 double]
    #     qs: [234x8 double]
    #        version: 'v 1.0.20.52'
    #        created: 'Wed Sep 02 11:48:41 202'
    #       DATE_OBS: '2015/12/18'
    #       TIME_OBS: '08:26:24.100'
    #         CDELT1: 2.8812163
    #        AZIMUTH: 10
    #        SOL_DEC: -23.379
    #        SOLAR_R: 975.25
    #        SOLAR_P: -351.10001
    #        SOLAR_B: -1.3
    #        RATAN_P: 13.271707
    #        N_FREQS: 8
    #       N_POINTS: 234
    # user_sel_freqs: '3.0000000e+000  4.0000000e+000  6.0000000e+000  8.0000000e+000  1.0000000e+001  1.2000000e+001  1.4000000e+001  1.6500000e+00'
    #          shift: 0


