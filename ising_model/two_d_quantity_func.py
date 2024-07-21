'''
동일한 온도, 격자 크기, 전이횟수에 대해서 알고리즘을 여러 번 실행하고, 그에 따른 변동을 확인하기 위함.
2차원 이징모형.
'''

import matplotlib.pyplot as plt
import numpy as np
import random
import datetime, time

np.set_printoptions(linewidth=np.inf)  # np.array의 모든 내용을 한 줄으로만 출력한다. 

class IsingModelFunc():
    def __init__(self):
        self.kb = 1
        self.J = 1
        self.L = 10
        self.T = 1
        self.transition = 10000
        self.lattice = []
        self.row, self.col = 0, 0

        self.all_eng = []
        self.all_mag = []

        self.init_latt = self.initialize_state()
        self.run = 10
        for r in range(self.run):
            self.lattice = self.init_latt.copy()  # 같은 초기 상태에 대해서 실험.
            energies = []
            magnets = []
            for tr in range(self.transition):
                self.metropolis_algorithm(self.T)
                energies.append(self.calc_current_energy())
                magnets.append(self.calc_current_mag())
            energy_mean = sum(energies) / (self.transition * (self.L**2))
            magnet_mean = sum(magnets) / (self.transition * (self.L**2))
            self.all_eng.append(energy_mean)
            self.all_mag.append(magnet_mean)

    def initialize_state(self):
        spin = [-1, 1]
        lattice = np.array([[random.choice(spin) for i in range(self.L)] for j in range(self.L)], dtype=np.float64)
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

    def metropolis_algorithm(self, t):
        self.pick_spin()
        self.calc_energy_diff()
        self.flip_spin(t)

    def calc_current_energy(self):
        curr_energy = 0
        row_energy = 0
        col_energy = 0
        for i in range(self.L):
            for j in range(self.L):
                row_energy = self.lattice[i][j] * self.lattice[(i+1) % self.L][j]
                col_energy = self.lattice[i][j] * self.lattice[i][(j+1) % self.L]
                curr_energy += row_energy + col_energy
        return (-self.J) * curr_energy

    def calc_current_mag(self):
        # return abs(np.sum(self.lattice))
        return abs(self.lattice.sum())

    def get_lattice(self):
        return self.lattice

    def show_quantity_graph(self):
        run_array = np.arange(1, self.run+1, 1)
        fig = plt.figure()
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        ax1 = fig.add_subplot(121)
        ax1.set_title('mean energy. size={}, tr={}, T={}'.format(self.L**2, self.transition, self.T))
        #ax1.set_ylim(-1, 0)
        ax1.set_xlabel('Try')
        ax1.set_ylabel('Energy per spin')
        ax1.plot(run_array, self.all_eng)
        ax2 = fig.add_subplot(122)
        ax2.set_title('mean magnetization')
        ax2.set_xlabel('Try')
        ax2.set_ylabel('Magnetization per spin')
        #ax2.set_ylim(0, 1)
        ax2.plot(run_array, self.all_mag)
        plt.show()


ising = IsingModelFunc()
ising.show_quantity_graph()