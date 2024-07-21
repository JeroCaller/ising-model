'''
1차원 이징 모형에 대해서, 
시각화 따로, 물리량 구하기 따로 했던 기존의 방식에서
magnetic domain 시각화와 물리량의 계산을 동시에 수행하는 방법으로 변경. 
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
        self.L = 1000
        self.T = np.arange(0.1, 4.0 + 0.1, 0.1)
        self.transition = 100000
        self.lattice = []
        self.col = 0

        self.all_eng = []
        self.all_mag = []
        sample_eng = []  # 각 전이 상태에서의 에너지. 각 온도마다 구분한다. 
        sample_mag = []  # 각 전이 상태에서의 자기화

        #self.sample_show_T = [0.1, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5] # 여러 온도 중 시각화로 보이고 싶은 온도를 표시한다.
        self.sample_show_T = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.1]
        self.sample_latt = []  # sample T에 대응되는 격자 추출.

        self.all_latt = dict() # 온도 별 격자 배열 상태 모두 저장. 온도: 격자상태로 표현. 

        start_time = time.time()  # 알고리즘 수행 시간 측정용

        self.init_latt = self.initialize_state()
        for t in self.T:
            self.lattice = self.init_latt.copy()
            energies = []
            magnets = []
            for tr in range(self.transition):
                self.metropolis_algorithm(t)
                energies.append(self.calc_current_energy())
                magnets.append(self.calc_current_mag())
            energy_mean = sum(energies) / (self.transition * (self.L))
            magnet_mean = sum(magnets) / (self.transition * (self.L))
            self.all_eng.append(energy_mean)
            self.all_mag.append(magnet_mean)
            sample_eng.append(energies)
            sample_mag.append(magnets)
            self.all_latt['{:.1f}'.format(t)] = self.lattice.copy()

        end_time = time.time()
        sec = end_time - start_time
        result_list = str(datetime.timedelta(seconds=sec)).split(".")
        print(result_list[0])
        self.algorithm_time = result_list[0]  # 알고리즘 수행 시간 측정용

        self.eng_std_err = self.calc_all_std_error(sample_eng) 
        self.mag_std_err = self.calc_all_std_error(sample_mag)

        self.extract_sample_latt()

    def initialize_state(self):
        '''
        각 스핀마다 spin up 상태일 확률 90%, spin down 상태일 확률 10% 부여. 
        '''
        lattice = []
        for i in range(self.L):
            rand = random.random()
            if rand <= 0.9:
                spin = 1
            else:
                spin = -1
            lattice.append(spin)
        print(np.array(lattice))
        return np.array(lattice)

    def initialize_state2(self):
        '''
        초기 스핀 정렬 방향을 완전 랜덤으로 설정한다. 
        '''
        spin = [-1, 1]
        lattice = np.array([random.choice(spin) for i in range(self.L)], dtype=np.float64)
        return lattice

    def pick_spin(self):
        col = random.randint(0, self.L-1)
        self.col = col

    def calc_energy_diff(self):
        curr_spin = self.lattice[self.col]
        # 주기 경계 조건을 도입한다.
        left = self.lattice[(self.col-1) % self.L]
        right = self.lattice[(self.col+1) % self.L]
        
        init_energy = -self.J * (left + right) * curr_spin
        fin_energy = -self.J * (left + right) * (-curr_spin)
        self.ediff =  fin_energy - init_energy

    def flip_spin(self, T_):
        if self.ediff <= 0:
            # 에너지 변화량이 음수거나 0일 경우 무조건 해당 스핀을 뒤집는다. 
            self.lattice[self.col] = -self.lattice[self.col]
        else:
            rand = random.random()  # 0, 1사이의 실수를 무작위로 추출.
            if rand <= np.exp(-self.ediff/(self.kb * T_)):
                self.lattice[self.col] = -self.lattice[self.col]
            else:
                return  # 스핀을 뒤집지 않는다.

    def metropolis_algorithm(self, t):
        self.pick_spin()
        self.calc_energy_diff()
        self.flip_spin(t)

    def calc_current_energy(self):
        curr_energy = 0
        col_energy = 0
        for j in range(self.L):
            col_energy = self.lattice[j] * self.lattice[(j+1) % self.L]
            curr_energy += col_energy
        return (-self.J) * curr_energy

    def calc_current_mag(self):
        # return abs(np.sum(self.lattice))
        return abs(self.lattice.sum())

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
        ax1.set_title('mean energy. size={}, transition={}'.format(self.L, self.transition))
        ax1.set_ylim(-1, 0)
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

        size = str(self.L)
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

    def show_ising_model(self):
        '''
        magnetic domain을 시각화
        '''
        fig = plt.figure()
        # fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.35)
        # left, bottom, right, top: 서브플롯 4면의 위치 조정. wspace, hspace 포함 모두 0에서 1사이의 값만 가능한 걸로 보임. 
        fig.subplots_adjust(wspace=0.1, hspace=0.3)  # wspace: 서브플롯 간 너비 간격 조절. hspace: 서브플롯 간 높이 조절. 
        cmap = plt.get_cmap('binary')
        fig_row = 2
        fig_col = 4
        ax_init = fig.add_subplot(fig_row, fig_col, 1)
        ax_init.set_title('Initial state, iteration={}, size={}'.format(self.transition, self.L))
        ax_init.get_yaxis().set_visible(False)
        ax_init.set_xlim(1, self.L)
        ax_init.matshow(self.init_latt.reshape(1, -1), cmap=cmap, aspect='auto', vmin=-1, vmax=1)
        ax = []
        for i in range(len(self.sample_show_T)):
            ax.append(fig.add_subplot(fig_row, fig_col, i+2))
            ax[i].set_title('T={} [K]'.format(self.sample_show_T[i]))
            ax[i].get_yaxis().set_visible(False)
            ax[i].set_xlim(1, self.L)
            ax[i].matshow(self.sample_latt[i].reshape(1, -1), cmap=cmap, aspect='auto', vmin=-1, vmax=1)
        plt.show()

    def get_lattice(self):
        return self.lattice

    def extract_sample_latt(self):
        '''
        모든 온도 구간 중 보고 싶은 온도의 격자만 추출한다. 
        '''
        #print(self.all_latt)
        for t_value in self.sample_show_T:
            self.sample_latt.append(self.all_latt[str(t_value)])
        #print(self.sample_latt)


ising = IsingModel()
ising.show_quantity_graph()
ising.save_data('one_d_quantity.txt', 'size,transition,cost_time,T,energy_m,magnet_m,std_err_eng,std_err_mag')
ising.show_ising_model()
print(ising.get_lattice())