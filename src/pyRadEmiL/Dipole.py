import numpy as np
import sunpy.sun.constants as sun
import astropy.units as u

class Dipole:
    POS_CENTER = 1
    POS_PHTSPH = 2
    #-------------------------------------------------------------------------------
    def __init__(self, field_ph = 3000, depth = 16
               , dipole_dir = [0, 1]
               , pos = POS_PHTSPH
                ):
        self.__field_ph = field_ph
        self.__depth = (depth*u.Mm).to_value(u.Rsun) if pos == self.POS_PHTSPH else 1
        self.__mu = 0.5*field_ph*depth**3
        self.__dipole_dir = dipole_dir

        pass

    def get_field(self, x, y, z):
        rxy = np.sqrt(x**2 + y**2);
        d = self.__depth if pos == POS_PHTSPH else 0
        r = np.sqrt(rxy**2 + d**2);
        cost = d/r;
        Bz = self.__mu/r^3 * (3*cost^2 - 1);
        if rxy != 0:
            sint = rxy/r;
            Btr = self.__mu/r^3 * (3*cost*sint);
            Bx = Btr * x/rxy;
            By = Btr * y/rxy;
        else:
            Bx = 0;
            By = 0;

        return dict(Bx = Bx, By = By, Bz = Bz)
