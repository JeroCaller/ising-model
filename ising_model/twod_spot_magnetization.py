import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['font.family'] = 'Malgun Gothic'  # 한글 글꼴 적용

fig, ax = plt.subplots()
ax.spines[['top', 'right']].set_visible(False)
plt.xlabel('T')
plt.ylabel('m (h=0)', labelpad=20)
# plt.title('2차원 이징모형의 자발적 자기화 m')
xlim_num, ylim_num = 3, 2
ax.set_xlim(0, xlim_num)
ax.set_ylim(0, ylim_num)

x_T = np.linspace(0.01, 100, 1000000)
m = np.power(1-(np.power(1 - np.square(np.tanh(1/x_T)), 4) / (16 * np.power(np.tanh(1/x_T), 4))), 1/8)

ms = m[np.isnan(m) == False]           # 배열 값이 none(nan)인 요소는 제외하고 숫자가 존재하는 배열까지만 추출
xy_arr = np.array(list(zip(x_T, ms)))  # x에 대응하는 y값을 (x, y) 형태로 짝지음.
m_min_idx = np.where(xy_arr[:, 1] == ms.min(axis=0))  # y(m)의 최소값의 xy_arr에서의 인덱스 찾기
xy_min = xy_arr[m_min_idx][0]                         # y(m)의 최소값에 대응하는 x값과 짝 지음.

ax.plot(x_T, m)
ax.plot(xy_min[0], xy_min[1], marker='o')   # xy_min[0] = x, xy_min[1] = y
plt.vlines(xy_min[0], 0, xy_min[1], color='gray', linestyle='dashed')
ax.text(xy_min[0] + 0.1 , xy_min[1], '({:.2f}, {:.2f})'.format(xy_min[0], xy_min[1]))
plt.show()