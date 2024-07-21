import matplotlib.pyplot as plt
import numpy as np

temp = np.arange(0.1, 3 + 0.1, 0.1)
S = np.log(2 * np.cosh(1/temp)) - np.tanh(1/temp) / temp

plt.plot(temp, S)
plt.xlabel('T [K]')
plt.ylabel('S/N')
plt.show()