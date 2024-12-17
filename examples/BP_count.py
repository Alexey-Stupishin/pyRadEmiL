def gauss(x, sigma):
    return np.exp(-((x / (sigma / (2 * np.sqrt(2 * np.log(2))))) ** 2) / 2) / (
        sigma / (2 * np.sqrt(2 * np.log(2)))  * np.sqrt(2 * np.pi) 
    )


import scipy

si = np.array([
    252,
    209.82,
    180,
    157.37,
    140,
    125.9,
    114.45,
    104.91,
    96.84,
    90,
    84,
    71.94,
    63.00,
    56,
    50.36,
    45.8,
    42,
    38.74,
    36.0,
    31.47,
    28,
    25.2,
    22.89,
    20.98,
    19.37,
    18.0,
    16.78,
    15.737,
    14.81,
    14.00,
]) / arcsecpersample 

measured_frqs = np.concatenate((np.arange(1, 3.2, .2), np.arange(3.5, 7.5, .5), np.arange(8, 19, 1))) 


def sigma(f):

    return np.interp(f, measured_frqs, si)
import numpy as np
import plotly.graph_objects as go

sun_size = 4095      
arcsecpersample = 0.5


import scipy.constants

INTEGRATION_INTERVAL_WIDTH = 60000 
APERTURE_WIDTH = 300 
ILLUMINATION_SCALING = 800
NPOINTS = 4095

x = np.arange(0,  sun_size)

def rectangle(x):
    return np.where(abs(x)<=0.5, 1, 0)

wl = np.zeros((num_frqs))

for i, e in enumerate(freqs_sc):
    wl[i] = scipy.constants.speed_of_light / freqs_sc[i] / 1e9



def twosinc1(x, freq):
    xx = x / 2 / 3600 / 180 * np.pi
    k = 2 * np.pi * freq * 1e9 / (scipy.constants.speed_of_light)
    y = k * xx
    a = APERTURE_WIDTH
    b = a / (2 * np.arccos(.2))
    return freq * (2*b*(-np.cos((a*y)/2)*np.sin(a/(2*b)) + b*y*np.cos(a/(2*b))*np.sin((a*y)/2)))/(-1 + b**2*y**2)


def twosinc_sq1(x, e):                                  
    return np.abs(twosinc1(x, e)) ** 2/e  / np.sum(np.abs(twosinc1(x, e)) ** 2/e )
    
    
dlr = []
drr = []
scalel = []
scaler = []


sc = .000025

for i, e in enumerate(freqs_sc):
    y = twosinc_sq1(x - sun_size//2 , e)
    # y = gauss(x - sun_size // 2, sigma(e))
    # y = y / y.max()

    dlr.append(np.convolve(y, sdls1[i] + sdls2[i] + sdls3[i], mode="same") )
    scalel.append(sc * dlr[-1].mean())
    noise = np.random.normal(loc=0, scale=scalel[-1], size=len(dlr[-1]))
    dlr[-1] = dlr[-1] + noise
    
    drr.append(np.convolve(y, sdrs1[i] + sdrs2[i] + sdrs3[i], mode="same") )
    scaler.append(sc * drr[-1].mean())
    noise = np.random.normal(loc=0, scale=scaler[-1], size=len(drr[-1]))
    drr[-1] += noise


dlrr = np.array(dlr)
drrr = np.array(drr)
