'''
2차원 이징 모형에 대해서, 
시각화 따로, 물리량 구하기 따로 했던 기존의 방식에서
magnetic domain 시각화와 물리량의 계산을 동시에 수행하는 방법으로 변경. 

ver2
- 연산 시간을 조금이라도 빠르게 하기 위한 최적화 진행 중. 
- 에러바 도입.
- capacity, magnet susceptibility(자화율) 물리량 추가 도입.

수정해야할 점
- 열용량이 온도 1도 이하 부근에서 부자연스럽게 혼자 높은 값을 차지함. 이유 파악 중. 
'''

from re import I
import matplotlib.pyplot as plt
import numpy as np
import random
import datetime, time

np.set_printoptions(linewidth=np.inf)  # np.array의 모든 내용을 한 줄으로만 출력한다. 
np.set_printoptions(precision=12, suppress=True) # 지수의 과학적표기법 대신 정확한 실수 표현으로 바꾼다. precision = 소수점 ?자리까지 나타낸다.

class IsingModel():
    def __init__(self):
        self.kb = 1
        self.J = 1
        self.L = 2
        #self.T = np.arange(0.1, 4.0 + 0.1, 0.1)
        #self.T = np.arange(1.00, 4.00 + 0.01, 0.01)
        self.T = np.round(np.linspace(0.1, 4.0, 50), 2)
        self.transition = (self.L * self.L) * 100  # 추출 수는 비용시간 및 오차 고려해서 적절히 선택.
        self.lattice = []
        self.row, self.col = 0, 0

        self.all_eng = []
        self.all_mag = []
        self.all_eng_sqr_in = [] # <o^2>
        self.all_mag_sqr_in = []
        self.all_eng_sqr_out = []  # <o>^2
        self.all_mag_sqr_out = []
        self.all_capacity = [] # 모든 온도에 대한 열용량
        self.all_suscep = [] # 모든 온도에 대한 자화율
        sample_eng = []  # 각 전이 상태에서의 에너지. 각 온도마다 구분한다. 
        sample_mag = []  # 각 전이 상태에서의 자기화

        self.sample_show_T = [0.1, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5] # 여러 온도 중 시각화로 보이고 싶은 온도를 표시한다. 7개.
        #self.sample_show_T = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.1]
        self.sample_latt = []  # sample T에 대응되는 격자 추출.

        self.all_latt = dict() # 온도 별 격자 배열 상태 모두 저장. 온도: 격자상태로 표현. 

        start_time = time.time()  # 알고리즘 수행 시간 측정용

        self.init_latt = self.initialize_state()
        for t in self.T:
            self.lattice = self.init_latt.copy()
            energies = []
            energies_sqr_in = []
            magnets = []
            magnets_sqr_in = []
            for tr in range(self.transition):
                self.metropolis_algorithm(t)
                current_energy = self.calc_current_energy()
                energies.append(current_energy)
                energies_sqr_in.append(current_energy ** 2)
                current_magnet = self.calc_current_mag()
                magnets.append(current_magnet)
                magnets_sqr_in.append(current_magnet ** 2)

            energy_mean = np.mean(energies)
            magnet_mean = np.mean(magnets)
            energy_mean_sqr_out = energy_mean ** 2 # <o>^2
            magnet_mean_sqr_out = magnet_mean ** 2
            energy_mean_sqr_in = np.mean(energies_sqr_in)
            magnet_mean_sqr_in = np.mean(magnets_sqr_in)
            self.all_eng.append(energy_mean / (self.L**2))
            self.all_mag.append(magnet_mean / (self.L**2))
            self.all_eng_sqr_in.append(energy_mean_sqr_in)
            self.all_mag_sqr_in.append(magnet_mean_sqr_in)
            self.all_eng_sqr_out.append(energy_mean_sqr_out)
            self.all_mag_sqr_out.append(magnet_mean_sqr_out)
            self.all_capacity.append(self.calc_capacity(t, energy_mean_sqr_in, energy_mean_sqr_out))
            self.all_suscep.append(self.calc_suscep(t, magnet_mean_sqr_in, magnet_mean_sqr_out))
            sample_eng.append(energies)
            sample_mag.append(magnets)
            self.all_latt['{:.1f}'.format(t)] = self.lattice.copy()

        end_time = time.time()
        time_interval = end_time - start_time
        cost_time_list = str(datetime.timedelta(seconds=time_interval)).split(".")
        print(cost_time_list[0])
        self.algorithm_time = cost_time_list[0]  # 알고리즘 수행 시간 측정용

        self.eng_std_err = self.calc_all_std_error(sample_eng) 
        self.mag_std_err = self.calc_all_std_error(sample_mag)
        self.capacity_std_err = self.calc_all_std_err_cx(self.all_eng_sqr_in, self.all_eng_sqr_out, 'c')
        self.suscep_std_err = self.calc_all_std_err_cx(self.all_mag_sqr_in, self.all_mag_sqr_out, 'x')

        self.tc_capacity = self.extract_tc_quantity(self.all_capacity)
        self.tc_suscep = self.extract_tc_quantity(self.all_suscep)

        self.extract_sample_latt()

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
        print(np.array(lattice))
        return np.array(lattice)

    def initialize_state2(self):
        '''
        초기 스핀 정렬 방향을 완전 랜덤으로 설정한다. 
        '''
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

        spin_sum = top + bottom + left + right
        self.ediff = 2 * self.J * spin_sum * curr_spin

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
        return abs(self.lattice.sum())

    def calc_capacity(self, t, sqr_in, sqr_out):
        # 한 온도에 대해서만 연산
        # t = 특정 온도
        c_one_temp = (sqr_in - sqr_out) / ((t**2) * (self.L**2))
        return c_one_temp

    def calc_suscep(self, t, sqr_in, sqr_out):
        # 한 온도에 대해서만 연산
        x_one_temp = (sqr_in - sqr_out) / (t * (self.L**2))
        return x_one_temp

    def extract_tc_quantity(self, quantity_list=[]):
        '''
        모든 온도에 대한 capacity 또는 suscep 그래프에서, 가장 높은 capacity 값을 가지는 임계온도의 수치값 추출
        capacity의 최대값에 해당하는 온도(임계온도)를 구하는 메소드.
        quantity = capacity or suscep (self.all_capacity or self.all_suscep)
        '''
        temps_list = list(np.round(self.T.copy(), 2))
        for i, value in enumerate(temps_list):
            if value >= 1.0:
                index_one_t = i
                break
        q_list = quantity_list.copy()
        temps_list = temps_list[index_one_t:]
        q_list = q_list[index_one_t:]
        max_quantity = max(q_list)
        max_index = q_list.index(max_quantity)
        Tc = temps_list[max_index]
        return Tc

    def calc_all_std_err_cx(self, q_sqr_in, q_sqr_out, mode):
        '''
        열용량 또는 자화율의 표준오차 구하는 함수.
        q_sqr_in or out = 열용량 또는 자화율의 제곱의 기대값 또는 기대값의 제곱 리스트 (모든 온도에 대해)
        mode)
        'c': capacity
        'x': magnetic susceptibility
        '''
        if len(q_sqr_in) != len(q_sqr_out):
            print('에러: 두 리스트의 길이가 맞지 않음.')
            return
        pair_q = list(zip(q_sqr_in, q_sqr_out))
        std_err = []
        for i, values in enumerate(pair_q):
            if mode == 'c':
                std_err.append(np.sqrt(values[0] - values[1]) / (np.sqrt(self.transition) * (self.T[i] ** 2)) )
            else:
                # 'x'
                std_err.append(np.sqrt(values[0] - values[1]) / ( np.sqrt(self.transition) * (self.T[i]) ) )
        return std_err

    def calc_all_std_error(self, all_samples):
        '''
        한 온도에 대한 표준 오차를 모든 온도에 대해 구한다. 
        에너지 및 자화 전용
        '''
        all_std_error = []
        for sample in all_samples:
            calc_error = np.std(sample) / np.sqrt(self.transition)
            all_std_error.append(calc_error)
        return all_std_error

    def show_quantity_graph(self, mode=''):
        fig = plt.figure()
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        ax1 = fig.add_subplot(221)
        ax1.set_title('mean energy. size={}, transition={}'.format(self.L**2, self.transition))
        ax1.set_ylim(-2, 0)
        ax1.set_xlabel('T [K]')
        ax1.set_ylabel('<E>/N')
        #ax1.plot(self.T, self.all_eng)
        #ax1.errorbar(self.T, self.all_eng, yerr=self.eng_std_err, capsize=3)

        ax2 = fig.add_subplot(222)
        ax2.set_title('mean magnetization')
        ax2.set_xlabel('T [K]')
        ax2.set_ylabel('<m>/N')
        ax2.set_ylim(0, 1)
        #ax2.plot(self.T, self.all_mag)
        #ax2.errorbar(self.T, self.all_mag, yerr=self.mag_std_err, capsize=3)

        ax3 = fig.add_subplot(223)
        ax3.set_title('specific heat capacity, Tc={} [K]'.format(round(self.tc_capacity, 2)))
        ax3.set_xlabel('T [K]')
        ax3.set_ylabel('C/N')
        #ax3.plot(self.T, self.all_capacity)
        #ax3.errorbar(self.T, self.all_capacity, yerr=self.capacity_std_err, capsize=3)

        ax4 = fig.add_subplot(224)
        ax4.set_title('magnetic susceptibility Tc={} [K]'.format(round(self.tc_suscep, 2)))
        ax4.set_xlabel('T [K]')
        ax4.set_ylabel('x/N')
        #ax4.plot(self.T, self.all_suscep)
        #ax4.errorbar(self.T, self.all_suscep, yerr=self.suscep_std_err, capsize=3)

        if mode == 'plot':
            ax1.plot(self.T, self.all_eng, 'o--')
            ax2.plot(self.T, self.all_mag, 'o--')
            ax3.plot(self.T, self.all_capacity, 'o--')
            ax4.plot(self.T, self.all_suscep, 'o--')
        else:
            # 'errorbar'
            ax1.errorbar(self.T, self.all_eng, yerr=self.eng_std_err, capsize=3, fmt='.--')
            ax2.errorbar(self.T, self.all_mag, yerr=self.mag_std_err, capsize=3, fmt='.--')
            ax3.errorbar(self.T, self.all_capacity, yerr=self.capacity_std_err, capsize=3, fmt='.--')
            ax4.errorbar(self.T, self.all_suscep, yerr=self.suscep_std_err, capsize=3, fmt='.--')        
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
        capacity = str(np.array(self.all_capacity))
        suscep = str(np.array(self.all_suscep))
        std_err_cap = str(np.array(self.capacity_std_err))
        std_err_sus = str(np.array(self.suscep_std_err))
        tc_cap = str(self.tc_capacity)
        tc_sus = str(self.tc_suscep)
        cost_time = str(self.algorithm_time)
        with open(file_name, 'a') as f:
            f.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(size, transition, cost_time, tc_cap, tc_sus, T, energy_m, magnet_m, std_err_eng, std_err_mag, capacity, suscep,\
                std_err_cap, std_err_sus))
        print('Data is successfully saved!')

    def save_all_latt(self, file_name, field):
        '''
        격자 배열 상태 데이터 저장
        '''
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

        init_latt = self.init_latt.reshape(1, -1)
        all_latt = dict()
        for key, value in self.all_latt.items():
            all_latt[key] = str(value.reshape(1, -1))

        with open(file_name, 'a') as f:
            f.write('{}@{}@{}\n'.format(self.T, init_latt, all_latt))
        print('lattice info is successfully saved!')

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
        ax_init.set_title('Initial state, iteration={}'.format(self.transition))
        ax_init.matshow(self.init_latt, cmap=cmap, vmin=-1, vmax=1)  # 각각 히트맵에 표현할 최소값, 최대값. 모든 값이 같은 값일 때 설정 안하면 오류남.
        ax = []
        for i in range(len(self.sample_show_T)):
            ax.append(fig.add_subplot(fig_row, fig_col, i+2))
            ax[i].set_title('T={} [K]'.format(self.sample_show_T[i]))
            ax[i].matshow(self.sample_latt[i], cmap=cmap, vmin=-1, vmax=1)
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


ising = IsingModel()
ising.show_quantity_graph('plot')
ising.show_quantity_graph('errorbar')
ising.save_data('two_d_quantity.txt', 'size,transition,cost_time,tc_cap,tc_sus,T,energy_m,magnet_m,std_err_eng,std_err_mag,capacity,suscep,std_err_cap,std_err_sus')
ising.save_all_latt('two_d_lattice.txt', 'T,init_latt,lattice')  # 현재 제대로 작동 안되는 상황
ising.show_ising_model()
#print(ising.get_lattice())