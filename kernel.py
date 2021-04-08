from math import pi
from scipy import integrate


class WendlandC6:

    def __init__(self, dimension=3):
        self.support_radius = 2.5
        self.dimension = dimension

    def __reformat(self, dr, h):
        dr = abs(dr)  # dr is usually a 3D vector, which is normed by dr = sqrt(r^2)
        H = self.support_radius * h
        s = dr / H
        return dr, H, s

    def W(self, dr, h):
        dr, H, s = self.__reformat(dr=dr, h=h)
        coeff = (1365.0 / 64.0) * (1.0 / (pi * (H ** self.dimension)))
        coeff2 = (1 - s) ** 8
        if (1 - s) < 0:
            coeff2 = 0
        w = coeff * (coeff2 * (1.0 + (8.0 * s) + (25.0 * s ** 2) + (32.0 * s ** 3)))
        return w

    def intW(self, dr, h):
        dr, H, s = self.__reformat(dr=dr, h=h)
        coeff = (1 / (64 * pi * (H ** self.dimension))) * (91 * s)
        intw = coeff * ((40 * s ** 11) - (315 * s ** 10) + (1056 * s ** 9) - (1925 * s ** 8) + (1980 * s ** 7) - (
                990 * s ** 6) + (198 * s ** 4) - (55 * s ** 2) + 15)
        return intw

    def __numericalW(self, s, dr, h):
        dr, H, s_redun = self.__reformat(dr=dr, h=h)
        coeff = (1365.0 / 64.0) * (1.0 / (pi * (H ** self.dimension)))
        coeff2 = (1 - s) ** 8
        if (1 - s) < 0:
            coeff2 = 0
        w = coeff * (coeff2 * (1.0 + (8.0 * s) + (25.0 * s ** 2) + (32.0 * s ** 3)))
        return w

    def numerical_intW(self, lower_bound, dr, h):
        dr, H, s = self.__reformat(dr=dr, h=h)
        upper_bound = s
        return integrate.quad(self.__numericalW, lower_bound, upper_bound, args=(dr, h,))[0]

    def mass_intW(self, dr, h):
        dr, H, s = self.__reformat(dr=dr, h=h)
        coeff = (1.0 / 16.0) * s ** 3
        return coeff * ((3120 * s ** 11) - (24255 * s ** 10) + (80080 * s ** 9) - (143325 * s ** 8)
                        + (144144 * s ** 7) - (70070 * s ** 6) + (12870 * s ** 4) - (3003 * s ** 2) + 455)
