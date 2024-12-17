import numpy as np
import scipy
import matplotlib.pyplot as plt

APERTURE_WIDTH = 300 
sun_size = 4095      

x = (np.arange(0,  sun_size) - sun_size//2)

def twosinc1(x, freq):
    xx = x / 3600 / 180 * np.pi
    k = 2 * np.pi * freq * 1e9 / (scipy.constants.speed_of_light)
    y = k * xx
    a = APERTURE_WIDTH
    b = a / (2 * np.arccos(.2))
    return freq * (2*b*(-np.cos((a*y)/2)*np.sin(a/(2*b)) + b*y*np.cos(a/(2*b))*np.sin((a*y)/2)))/(-1 + b**2*y**2)


def twosinc_sq1(x, e):                                  
    # return np.abs(twosinc1(x, e)) ** 2/e  / np.sum(np.abs(twosinc1(x, e)) ** 2/e )
    sq = twosinc1(x, e) ** 2
    return sq  / np.sum(sq)
    
def gauss(x, sigma):
    return np.exp(-((x / (sigma / (2 * np.sqrt(2 * np.log(2))))) ** 2) / 2) / (
        sigma / (2 * np.sqrt(2 * np.log(2)))  * np.sqrt(2 * np.pi) 
    )

freq = 3
y = twosinc_sq1(x, freq)
y2 = gauss(x, 84)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x, y)
ax.plot(x, y2)
plt.xlim([-300, 300])
plt.show(block = False)

pass
