'''
메트로폴리스 알고리즘을 이용한 2차원 이징모형의 시각화. 
'''


import matplotlib.pyplot as plt
import numpy as np
import random

class IsingModel():
    def __init__(self):
        self.J = 1
        self.kb = 1
        self.T = 5
        self.L = 300
        self.iteration = (self.L ** 2) * 100
        self.lattice = []
        self.row, self.col = 0, 0
        self.ediff = 0

        self.initialize_state()
        self.init_latt = self.lattice.copy()
        for i in range(self.iteration):
            self.pick_spin()
            self.calc_energy_diff()
            self.flip_spin()

    def initialize_state(self):
        '''
        2차원 격자의 좌상칸을 초기점 0, 0으로 잡는다. 오른쪽으로 갈수록 열이, 아래로 갈수록 행의 인덱스가 증가하는 식이다. 
        '''
        spin = [-1, 1]
        lattice = np.array([[random.choice(spin) for i in range(self.L)] for j in range(self.L)])
        self.lattice = lattice

    def get_lattice(self):
        return self.lattice

    def pick_spin(self):
        row = random.randint(0, self.L-1)
        col = random.randint(0, self.L-1)
        self.row, self.col = row, col

    def calc_energy_diff(self):
        curr_spin = self.lattice[self.row][self.col]
        top, bottom, left, right = 0, 0, 0, 0
        # 주기 경계 조건을 도입한다.
        if self.row == 0:
            top = self.lattice[self.L-1][self.col]
        elif self.row == self.L-1:
            bottom = self.lattice[0][self.col]
        elif self.col == 0:
            left = self.lattice[self.row][self.L-1]
        elif self.col == self.L-1:
            right = self.lattice[self.row][0]
        else:
            top = self.lattice[self.row-1][self.col]
            bottom = self.lattice[self.row+1][self.col]
            left = self.lattice[self.row][self.col-1]
            right = self.lattice[self.row][self.col+1]
        
        init_energy = -self.J * (top + bottom + left + right) * curr_spin
        fin_energy = -self.J * (top + bottom + left + right) * (-curr_spin)
        self.ediff =  fin_energy - init_energy

    def flip_spin(self):
        if self.ediff <= 0:
            # 에너지 변화량이 음수거나 0일 경우 무조건 해당 스핀을 뒤집는다. 
            self.lattice[self.row][self.col] = -self.lattice[self.row][self.col]
        else:
            rand = random.random()  # 0, 1사이의 실수를 무작위로 추출.
            if rand <= np.exp(-self.ediff/(self.kb * self.T)):
                self.lattice[self.row][self.col] = -self.lattice[self.row][self.col]
            else:
                return  # 스핀을 뒤집지 않는다. 

    def show_ising_model(self):
        fig, ax = plt.subplots(1, 2)
        cmap = plt.get_cmap('binary')
        ax[0].set_title('initial state')
        ax[0].matshow(self.init_latt, cmap=cmap)
        ax[1].set_title('iteration={}, T={}'.format(self.iteration, self.T))
        ax[1].matshow(self.lattice, cmap=cmap)
        plt.show()

ising = IsingModel()
ising.show_ising_model()