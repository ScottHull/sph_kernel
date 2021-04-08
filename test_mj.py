from math import pi
import numpy as np
from kernel import WendlandC6
import matplotlib.pyplot as plt


def calcSmoothedMass(kernel, m_j, dr, h):
    # return m_j * kernel.mass_intW(dr=dr, h=h)
    H = kernel.support_radius * h
    return 4 * pi * (H ** 3) * m_j * kernel.intW_s2(dr=dr, h=h)


kernel = WendlandC6()

m_j = 10 ** 20
h = 1
dr_range = np.arange(-2.5, 2.5, 0.1)
masses = [calcSmoothedMass(kernel=kernel, m_j=m_j, dr=dr, h=h) for dr in dr_range]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(
    [abs(dr / h) for dr in dr_range],
    masses,
    linewidth=2.0,
    color='black'
)
ax.axhline(m_j, color='red', linewidth=2.0, linestyle="--", label="m_j")
ax.axvline(2 * h, color='blue', linewidth=2.0, linestyle="--", label="2h")
ax.set_xlabel("s = r / h")
ax.set_ylabel("m_j^*")
ax.set_title("Wendlend C6")
ax.grid()
ax.legend()
plt.show()
