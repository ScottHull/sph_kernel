import numpy as np
from kernel import WendlandC6
import matplotlib.pyplot as plt

kernel = WendlandC6()

dr_range = np.arange(-2.5, 2.5, 0.1)
# dr_range = [(1.62918e6, 565945), (1.36529e6, 565945)]
# expected values: 1.36529e-27
h = 1
calculated_w = [kernel.W(dr=dr, h=h) for dr in dr_range]
# calculated_w = [kernel.W(dr=dr[0], h=dr[1]) for dr in dr_range]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(
    [abs(dr / h) for dr in dr_range],
    # [dr[0] / dr[1] for dr in dr_range],
    calculated_w,
    linewidth=2.0,
    color='black'
)
ax.set_xlabel("s = r / h")
ax.set_ylabel("W(s)")
ax.set_title("Wendlend C6")
ax.grid()
plt.show()
