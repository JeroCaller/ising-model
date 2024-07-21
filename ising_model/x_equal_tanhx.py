'''
x = tanh(ax)의 그래프 그리기
'''

import matplotlib.pyplot as plt
import numpy as np


fig, ax = plt.subplots(1, 2)
for i in range(2):
    ax[i].spines[['left', 'bottom']].set_position(("data", 0))  # xy좌표축으로 변환하여 그림
    ax[i].spines[['top', 'right']].set_visible(False)
    ax[i].plot(1, 0, ">k", transform=ax[i].get_yaxis_transform(), clip_on=False)  # x축 화살표 그리기
    ax[i].plot(0, 1, "^k", transform=ax[i].get_xaxis_transform(), clip_on=False)  # y축 화살표 그리기
    ax[i].set_xlim(-2, 2)
    ax[i].set_ylim(-2, 2)
    ax[i].set_xlabel("m", loc='right')
    ax[i].set_ylabel('y', loc='top', rotation=0)

# plt.title(r'$m = \tanh (\beta nJm)$', loc='center', pad=20)
point_config = {'start': -10, 'end': 10, 'space': 10000}

# ====== 데이터 입력 ======
x1 = np.linspace(point_config['start'], point_config['end'], point_config['space'])
y1 = x1

rng = np.random.default_rng(1)
more_than_one = rng.choice(np.arange(2, 5))  # > 1
less_than_one = rng.random()

x2 = np.linspace(point_config['start'], point_config['end'], point_config['space'])
y2 = np.tanh(more_than_one * x2)

x3 = np.linspace(point_config['start'], point_config['end'], point_config['space'])
y3 = np.tanh(less_than_one * x2)
# ====== 데이터 입력 완료======
ax[0].plot(x1, y1, x2, y2)
ax[1].plot(x1, y1, x3, y3)
ax[0].text(-2, 2, r'$\beta qJ > 1$', fontdict={'size':15})
ax[1].text(-2, 2, r'$\beta qJ < 1$', fontdict={'size':15})
ax[0].text(1.4, 1.9, 'y = m')
ax[0].text(1.4, 0.8, r'$y = tanh(\beta qJm)$')
ax[1].text(1.4, 1.9, 'y = m')
ax[1].text(1.4, 0.7, r'$y = tanh(\beta qJm)$')
plt.show()