'''
2차원 이징 모형
초기 상태들을 무작위로 설정한 다음, 무작위로 임의의 칸을 선택한 후 뒤집는다.
'''

import matplotlib.pyplot as plt
import numpy as np
import random


def initialize_state(N=1):
    state = [-1, 1]
    lattice_data = []
    for i in range(N):
        one_row = [random.choice(state) for j in range(N)]
        lattice_data.append(one_row)
    return lattice_data


N = 100
lattice = np.array(initialize_state(N))
init_lattice = lattice.copy()

iteration = 10000
for iter in range(iteration):
    row = random.randint(0, N-1)
    col = random.randint(0, N-1)
    if lattice[row][col] == -1:
        lattice[row][col] = 1
    else:
        lattice[row][col] = -1

fig, ax = plt.subplots(1, 2)
cmap = plt.get_cmap('binary')  # 색깔을 흑백으로 바꿈. 
ax[0].set_title('initial lattice')
ax[0].matshow(init_lattice, cmap=cmap)
ax[1].set_title('result lattice')
ax[1].matshow(lattice, cmap=cmap)
plt.show()