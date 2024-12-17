import numpy as np
import plotly.graph_objects as go
import scipy
import scipy.constants

sun_size = 4095      
arcsecpersample = 0.5
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

INTEGRATION_INTERVAL_WIDTH = 60000 
APERTURE_WIDTH = 300 
ILLUMINATION_SCALING = 800
NPOINTS = 4095

x = np.arange(0,  sun_size)

measured_frqs = np.concatenate((np.arange(1, 3.2, .2), np.arange(3.5, 7.5, .5), np.arange(8, 19, 1))) 
num_frqs = len(measured_frqs)

print(measured_frqs)

pass
