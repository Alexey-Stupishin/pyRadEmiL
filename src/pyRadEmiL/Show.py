import matplotlib.pyplot as plt
import matplotlib.colors as clr
from matplotlib import cm
import numpy as np

class Show:
    @staticmethod
    def image(title, v, cblabel = None):
        plt.figure(title)
        c = plt.imshow(v, cmap ='Greens', interpolation ='nearest', origin ='lower')
        plt.colorbar(c, label = cblabel)
        plt.show(block = False)
        return c

    @staticmethod
    def image_log(title, v):
        plt.figure(title)
        c = plt.imshow(v, cmap ='Greens', interpolation ='nearest', origin ='lower', norm = clr.LogNorm(vmin = v.max()*1e-5, vmax = v.max()))
        plt.colorbar(c)
        plt.show(block = False)
        return c

    @staticmethod
    def plot(title, x, y, ax = None, color = 'k', linestyle = 'solid', linewidth = 2, label = None):
        if ax is None:
            fig = plt.figure(title)
            ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, y, color = color, linestyle = linestyle, linewidth = linewidth, label = label)
        plt.show(block = False)
        return ax

    @staticmethod
    def plot_y(title, y, ax = None):
        if ax is None:
            fig = plt.figure(title)
            ax = fig.add_subplot(1, 1, 1)
        ax.plot(range(len(y)), y)
        plt.show(block = False)
        return ax

    @staticmethod
    def surf(title, xg, yg, map, ax = None):
        fig, ax = plt.subplots(subplot_kw = {"projection": "3d"})
        surf = ax.plot_surface(xg, yg, map, cmap = cm.coolwarm,
                       linewidth=0, antialiased = False)
        
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()

        return ax
