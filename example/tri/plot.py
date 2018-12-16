#%matplotlib inline
#%config InlineBackend.figure_format = 'retina'

import numpy as np
import matplotlib.pyplot as plt

#%data01 = np.loadtxt("./result.txt", comments='#')

data01_axis1, data01_value1 = np.loadtxt("./data01.txt", comments='#', unpack=True)

fig = plt.figure(figsize=(4, 6))
ax = fig.add_subplot(111)
ax.plot(data01_axis1, data01_value1, "o-", color="k", label="value1 of data01")
ax.set_xlim(-0.5, 5.5)
ax.set_ylim(-0.5, 25.5)
ax.set_xlabel("axis1")
ax.set_ylabel("value1")
ax.legend(loc="upper left")
plt.show()
