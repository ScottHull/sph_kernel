import numpy as np
from kernel import WendlandC6
import matplotlib.pyplot as plt

kernel = WendlandC6()

dr_range = np.arange(-2.5, 2.5, 0.1)
h = 1
calculated_intw = [kernel.intW(dr=dr, h=h) for dr in dr_range]
numerical_intw = [kernel.numerical_intW(lower_bound=0, dr=dr, h=h) for dr in dr_range]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(
    [abs(dr / h) for dr in dr_range],
    calculated_intw,
    linewidth=2.0,
    color='black',
    label='analytical'
)
ax.plot(
    [abs(dr / h) for dr in dr_range],
    numerical_intw,
    linewidth=2.0,
    color='red',
    label='numerical',
    linestyle="--"
)
ax.set_xlabel("s = r / h")
ax.set_ylabel("int_0^s W(s)")
ax.set_title("Wendlend C6")
ax.grid()
ax.legend()
plt.show()
