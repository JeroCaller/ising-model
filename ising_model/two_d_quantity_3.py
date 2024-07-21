'''
ver3
오로지 전이횟수만 담겨져 있는 2차원 이징 모형의 물리량 측정 프로그램. 
전이횟수를 최대한으로 해서 물리량의 변동을 최대한 줄인다. 
'''

import matplotlib.pyplot as plt
import numpy as np
import random
import datetime, time

np.set_printoptions(linewidth=np.inf)  # np.array의 모든 내용을 한 줄으로만 출력한다. 

class IsingModel():
    def __init__(self):
        self.kb = 1
        self.J = 1
        self.L = 4
        self.T = np.arange(1.0, 4.0 + 0.1, 0.1)
        self.transition = 100000
        self.lattice = []
        self.row, self.col = 0, 0

        self.all_eng = []
        self.all_mag = []
        sample_eng = []  # 각 전이 상태에서의 에너지. 각 온도마다 구분한다. 
        sample_mag = []  # 각 전이 상태에서의 자기화

        start_time = time.time()

        self.init_latt = self.initialize_state()
        for t in self.T:
            self.lattice = self.init_latt.copy()
            energies = []
            magnets = []
            for tr in range(self.transition):
                self.metropolis_algorithm(t)
                energies.append(self.calc_current_energy())
                magnets.append(self.calc_current_mag())
            energy_mean = sum(energies) / (self.transition * (self.L**2))
            magnet_mean = sum(magnets) / (self.transition * (self.L**2))
            self.all_eng.append(energy_mean)
            self.all_mag.append(magnet_mean)
            sample_eng.append(energies)
            sample_mag.append(magnets)

        end_time = time.time()
        sec = end_time - start_time
        result_list = str(datetime.timedelta(seconds=sec)).split(".")
        print(result_list[0])
        self.algorithm_time = result_list[0]  # 알고리즘 수행 시간 측정용

        self.eng_std_err = self.calc_all_std_error(sample_eng) 
        self.mag_std_err = self.calc_all_std_error(sample_mag)

    def initialize_state(self):
        '''
        각 스핀마다 spin up 상태일 확률 90%, spin down 상태일 확률 10% 부여. 
        '''
        lattice = []
        for i in range(self.L):
            one_row_lattice = []
            for j in range(self.L):
                rand = random.random()
                if rand <= 0.1:
                    spin = 1
                else:
                    spin = -1
                one_row_lattice.append(spin)
            lattice.append(one_row_lattice)
        return np.array(lattice)

    def initialize_state2(self):
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

    def calc_std_error(self, sample):
        '''
        표본의 표준 편차(표준 오차)를 구한다. 한 온도에 대해서만 구하는 메소드.
        '''
        estimate_mean = sum(sample) / len(sample)  #물리량의 기대값과 동일
        variance = 0
        for one_sample in sample:
            variance += (one_sample-estimate_mean) ** 2
        variance = variance / (len(sample) * (len(sample) - 1))
        std_error = np.sqrt(variance)
        return std_error

    def calc_all_std_error(self, all_samples):
        '''
        한 온도에 대한 표준 오차를 모든 온도에 대해 구한다. 
        '''
        all_std_error = []
        for sample in all_samples:
            all_std_error.append(self.calc_std_error(sample))
        return all_std_error

    def show_quantity_graph(self):
        fig = plt.figure()
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        ax1 = fig.add_subplot(121)
        ax1.set_title('mean energy. size={}, transition={}'.format(self.L**2, self.transition))
        #ax1.set_ylim(-1, 0)
        ax1.set_xlabel('T [K]')
        ax1.set_ylabel('Energy per spin')
        ax1.plot(self.T, self.all_eng)
        ax2 = fig.add_subplot(122)
        ax2.set_title('mean magnetization')
        ax2.set_xlabel('T [K]')
        ax2.set_ylabel('Magnetization per spin')
        ax2.set_ylim(0, 1)
        ax2.plot(self.T, self.all_mag)
        plt.show()

    def save_data(self, file_name, field):
        try:
            with open(file_name, 'r') as f:
                data = f.readlines()
        except FileNotFoundError:
            with open(file_name, 'w') as f:
                f.write(field + '\n')  # 파일 생성
        else:
            if len(data) == 0 or data[0].strip() != field:
                with open(file_name, 'a') as f:
                    f.write(field + '\n')

        size = str(self.L ** 2)
        transition = str(self.transition)
        T = str(self.T)
        energy_m = str(np.array(self.all_eng))
        magnet_m = str(np.array(self.all_mag))
        std_err_eng = str(np.array(self.eng_std_err))
        std_err_mag = str(np.array(self.mag_std_err))
        cost_time = str(self.algorithm_time)
        with open(file_name, 'a') as f:
            f.write('{},{},{},{},{},{},{},{}\n'.format(size, transition, cost_time, T, energy_m, magnet_m, std_err_eng, std_err_mag))
        print('Data is successfully saved!')


ising = IsingModel()
ising.show_quantity_graph()
ising.save_data('two_d_quantity.txt', 'size,transition,cost_time,T,energy_m,magnet_m,std_err_eng,std_err_mag')