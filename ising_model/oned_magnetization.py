import matplotlib.pyplot as plt
import numpy as np
import random


fig, ax = plt.subplots()
ax.spines[['left', 'bottom']].set_position(("data", 0))  # xy좌표축으로 변환하여 그림
ax.spines[['top', 'right']].set_visible(False)
ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)  # x축 화살표 그리기
ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)  # y축 화살표 그리기
xlim_num, ylim_num = 10, 2
ax.set_xlim(-xlim_num, xlim_num)
ax.set_ylim(-ylim_num, ylim_num)

plt.xlabel('h', loc='right')
plt.ylabel('m', loc='top', rotation=0)

h = np.linspace(-10, 10, 100)
beta = 1 / random.randint(1, 10)
J = 1
m = np.sinh(beta * h) / np.sqrt(np.square(np.sinh(beta * h)) + np.exp(-4 * beta * J))

ax.plot(h, m)
plt.show()