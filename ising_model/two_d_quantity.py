import matplotlib.pyplot as plt
import numpy as np
import random
import datetime, time


class IsingModel():
    def __init__(self):
        self.kb = 1
        self.J = 1
        self.L = 16
        self.T = np.linspace(0.1, 5, 50)
        self.transition = 1000
        self.lattice = []
        self.init_latt = []
        self.row, self.col = 0, 0

        self.energy = 0
        self.magnet = 0
        self.all_eng = []
        self.all_mag = []

        self.init_latt = self.initialize_state()
        self.run = 10
        for t in self.T:
            energies = []
            magnets = []
            for r in range(self.run):
                self.lattice = self.init_latt.copy()
                self.energy = 0
                self.magnet = 0
                self.metropolis_algorithm(t)
                energies.append((-(self.J / self.transition) * self.energy) / (self.L**2))
                magnets.append((self.magnet / self.transition) / (self.L**2))
            self.all_eng.append(np.mean(energies))
            self.all_mag.append(np.mean(magnets))

    def initialize_state(self):
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

    def metropolis_algorithm(self, t):
        for i in range(self.transition):
            self.pick_spin()
            self.calc_energy_diff()
            self.flip_spin(t)
            self.energy += self.calc_current_energy()
            self.magnet += self.calc_current_mag()

    def calc_current_energy(self):
        curr_energy = 0
        row_energy = 0
        col_energy = 0
        for i in range(self.L):
            for j in range(self.L):
                row_energy = self.lattice[i][j] * self.lattice[(i+1) % self.L][j]
                col_energy = self.lattice[i][j] * self.lattice[i][(j+1) % self.L]
                curr_energy += row_energy + col_energy
        return curr_energy/2

    def calc_current_mag(self):
        # return abs(np.sum(self.lattice))
        return abs(self.lattice.sum())

    def show_quantity_graph(self):
        fig = plt.figure()
        fig.subplots_adjust(wspace=0.3, hspace=0.3)
        ax1 = fig.add_subplot(121)
        ax1.set_title('mean energy. size={}'.format(self.L**2))
        ax1.set_ylim(-1, 0)
        ax1.set_xlabel('T')
        ax1.set_ylabel('Energy per spin')
        ax1.plot(self.T, self.all_eng)
        ax2 = fig.add_subplot(122)
        ax2.set_title('mean magnetization')
        ax2.set_xlabel('T')
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
        run = str(self.run)
        T = str(self.T)
        energy_m = str(np.array(self.all_eng))
        magnet_m = str(np.array(self.all_mag))
        with open(file_name, 'a') as f:
            f.write('{},{},{},{},{},{}\n'.format(size, transition, run, T, energy_m, magnet_m))
        print('Data is successfully saved!')


start_time = time.time()

ising = IsingModel()
ising.show_quantity_graph()
#ising.save_data('two_d_quantity.txt', 'size,transition,run,T,energy_m,magnet_m')

end_time = time.time()
sec = end_time - start_time
result_list = str(datetime.timedelta(seconds=sec)).split(".")
print(result_list[0])