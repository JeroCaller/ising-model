'''
메트로폴리스 알고리즘을 이용한 2차원 이징모형의 시각화. 
온도의 변화에 따른 모습을 확인하기 위한 버전. 반복횟수, 격자 사이즈 고정. 
'''

import matplotlib.pyplot as plt
import numpy as np
import random

class IsingModel():
    def __init__(self):
        self.J = 1
        self.kb = 1
        #self.T = [1, 0.1, 0.01, 0.001, 0.0001] # 5개
        self.T = [1, 2, 2.27, 3, 5]
        self.L = 200
        self.iteration = (self.L ** 2) * 100
        self.lattice = []
        self.all_latt = []
        self.row, self.col = 0, 0
        self.ediff = 0

        self.init_latt = self.initialize_state()
        for t in self.T:
            self.lattice = self.init_latt.copy()
            for i in range(self.iteration):
                self.pick_spin()
                self.calc_energy_diff()
                self.flip_spin(t)
            self.all_latt.append(self.lattice)

    def initialize_state(self):
        '''
        각 스핀마다 spin up 상태일 확률 90%, spin down 상태일 확률 10% 부여. 
        '''
        lattice = []
        for i in range(self.L):
            one_row_lattice = []
            for j in range(self.L):
                rand = random.random()
                if rand <= 0.9:
                    spin = 1
                else:
                    spin = -1
                one_row_lattice.append(spin)
            lattice.append(one_row_lattice)
        return np.array(lattice)
        
    def initialize_state2(self):
        '''
        2차원 격자의 좌상칸을 초기점 0, 0으로 잡는다. 오른쪽으로 갈수록 열이, 아래로 갈수록 행의 인덱스가 증가하는 식이다. 
        '''
        spin = [-1, 1]
        lattice = np.array([[random.choice(spin) for i in range(self.L)] for j in range(self.L)])
        return lattice

    def pick_spin(self):
        row = random.randint(0, self.L-1)
        col = random.randint(0, self.L-1)
        self.row, self.col = row, col

    def calc_energy_diff(self):
        curr_spin = self.lattice[self.row][self.col]
        # 주기 경계 조건을 도입한다.
        top = self.lattice[(self.row-1) % self.L][self.col]
        bottom = self.lattice[(self.row+1) % self.L][self.col]
        left = self.lattice[self.row][(self.col-1) % self.L]
        right = self.lattice[self.row][(self.col+1) % self.L]
        
        init_energy = -self.J * (top + bottom + left + right) * curr_spin
        fin_energy = -self.J * (top + bottom + left + right) * (-curr_spin)
        self.ediff =  fin_energy - init_energy

    def flip_spin(self, T_):
        if self.ediff <= 0:
            # 에너지 변화량이 음수거나 0일 경우 무조건 해당 스핀을 뒤집는다. 
            self.lattice[self.row][self.col] = -self.lattice[self.row][self.col]
        else:
            rand = random.random()  # 0, 1사이의 실수를 무작위로 추출.
            if rand <= np.exp(-self.ediff/(self.kb * T_)):
                self.lattice[self.row][self.col] = -self.lattice[self.row][self.col]
            else:
                return  # 스핀을 뒤집지 않는다. 

    def show_ising_model(self):
        fig = plt.figure()
        # fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.35)
        # left, bottom, right, top: 서브플롯 4면의 위치 조정. wspace, hspace 포함 모두 0에서 1사이의 값만 가능한 걸로 보임. 
        fig.subplots_adjust(wspace=0.1, hspace=0.3)  # wspace: 서브플롯 간 너비 간격 조절. hspace: 서브플롯 간 높이 조절. 
        cmap = plt.get_cmap('binary')
        ax_init = fig.add_subplot(231)
        ax_init.set_title('Initial random state, iteration={}'.format(self.iteration))
        ax_init.matshow(self.init_latt, cmap=cmap)
        ax = []
        for i in range(len(self.T)):
            ax.append(fig.add_subplot(2, 3, i+2))
            ax[i].set_title('T={} [K]'.format(self.T[i]))
            ax[i].matshow(self.all_latt[i], cmap=cmap)
        plt.show()

    def get_lattice(self):
        return self.lattice


ising = IsingModel()
ising.show_ising_model()
print(ising.get_lattice())