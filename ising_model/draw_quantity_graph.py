'''
for 2d ising.
'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipk
from scipy.special import ellipe
import scipy.special
from scipy.optimize import fsolve

class ProcessData():
    def __init__(self, file_name, one_valued_field):
        self.file_name = file_name
        self.raw_data = []
        self.data = []
        self.field = []
        self.one_valued_field = one_valued_field # 필드 값 중 단일 값으로만 이루어진 필드의 개수

        self.read_data()
        self.process_data()
        self.extract_field()
        self.delete_field()
        self.process_data2()

    def read_data(self):
        with open(self.file_name, 'r') as f:
            self.raw_data = f.readlines()

    def process_data(self):
        self.data = [one_data.strip().split(',') for one_data in self.raw_data]

    def extract_field(self):
        self.field = {self.data[0][i]: i for i in range(len(self.data[0]))}

    def process_data2(self):
        for one_data in self.data:
            for i in range(self.one_valued_field, len(one_data)):
                one_data[i] = one_data[i][1:-2].split(' ')  # [, ] 제외
                for j in range(one_data[i].count('')):
                    one_data[i].remove('')
                for j, str_num in enumerate(one_data[i]):
                    one_data[i][j] = float(str_num)
                one_data[i] = np.array(one_data[i])

    def delete_field(self):
        del self.data[0]

    def get_data(self):
        return self.data

    def get_field(self):
        return self.field


class RearrangeData():
    def __init__(self, data, field, one_valued_data):
        self.data = data
        self.field = field
        self.one_valued_data = one_valued_data # 값이 리스트가 아닌 단일 값을 가지는 필드들의 개수
        
        self.rearr_data = []
        for i, one_data in enumerate(self.data):
            self.rearr_data.append(self.combine_data(one_data))  # 모든 격자에 대한 정보 포함.

    def zip_data(self, one_data):
        '''
        같은 격자 사이즈 내에서 각 온도에 대응되는 물리량 및 표준오차를 묶는다.
        '''
        data_set = list(zip(one_data[self.field['T']], one_data[self.field['energy_m']], one_data[self.field['magnet_m']], 
        one_data[self.field['std_err_eng']], one_data[self.field['std_err_mag']], one_data[self.field['capacity']], 
        one_data[self.field['suscep']], one_data[self.field['std_err_cap']], one_data[self.field['std_err_sus']]))  # 수정 필요!!!
        return data_set

    def combine_data(self, one_data):
        '''
        온도에 따른 물리량 및 표준오차로 묶인 데이터를 그 데이터에 해당하는 격자 사이즈, 전이횟수 등의 정보와 결합하여 하나의 데이터행 완성
        '''
        comb_data = [one_data[i] for i in range(self.one_valued_data)]
        comb_data.append(self.zip_data(one_data))
        return comb_data

    def get_re_data(self):
        return self.rearr_data

    def save_re_data(self, file_name):
        with open(file_name, 'w') as f:
            for i, one_data in enumerate(self.rearr_data):
                f.write('start{}\n'.format(i+1))
                for j in range(self.one_valued_data):
                    f.write('{},'.format(one_data[j]))
                f.write('\n')
                for one_set in one_data[-1]:
                    f.write('{}\n'.format(one_set))
                f.write('end{}\n'.format(i+1))

'''
def mean_field_magnet(q, T):
    Tc = q
    m = np.power(1-(T/Tc), 1/2)
    return m

def mean_field_energy(q, T):
    Tc = q
    energy_expectation = - Tc * np.power(np.tanh(np.sqrt((Tc/T)*(Tc/T-1))),2)
    return energy_expectation
'''
def exact_energy_2d(T):
    q = 2 * np.sinh(2/T) / np.power(np.cosh(2/T), 2)
    #print(q)
    #K1_approx = np.pi / 2 + (q**2)*(np.pi**3)/48 - (q**2) * (np.pi**5)/960 + (q**2)*(np.pi**7)/40320 - (q**2)*(np.pi**9)/2903040
    #K1_approx = np.pi / 2 + (q**2)*(np.pi ** 3) / 48 + (q**2) * (9 * (q ** 2) - 4) * (np.pi ** 5) / 3840
    #K1_approx = ellipk(q)
    K1_approx = scipy.special.ellipk(q**2)
    part1 = 1 / np.tanh(2/T)
    part2 = 2 * np.power(np.tanh(2/T), 2) - 1
    energy = - part1 * (1 + (2/np.pi) * part2 * K1_approx)
    return energy


class Newton_Lapson():
    def __init__(self, temps, nns):
        self.T = temps
        self.one_temp = 1
        self.nns = nns

        self.magnets = self.calc_many_solutions()
        self.energies = self.calc_mft_energy()

    def f(self, m):
        formula = np.tanh((self.nns/self.one_temp)*m) - m
        return formula

    def calc_many_solutions(self):
        solutions = []
        for t in self.T:
            self.one_temp = t
            solutions.append(fsolve(self.f, [1, 100])[0])
        return solutions

    def calc_mft_energy(self):
        energy = - (self.nns/2) * np.power(self.magnets, 2)
        return energy

    def get_mfa_magnets(self):
        return self.magnets

    def get_mfa_energies(self):
        return self.energies


class DrawGraph():
    def __init__(self, data, field, dim):
        self.data = data
        self.field = field
        self.temp = self.data[1][self.field['T']]
        self.dim = dim #1차원 또는 2차원. 차원 수 
        self.linestyles = ['o--', 'v--', '^--', 's--', 'D--', '*--', 'p--']

    def show_graph(self, plot_mode):
        # plot_mode => 'plot', 'errorbar' 두 개 가능
        if self.dim == 1:
            mfa_inst = Newton_Lapson(self.temp, 2)
            mfa_magnets = mfa_inst.get_mfa_magnets()
            mfa_energies = mfa_inst.get_mfa_energies()
            exact_energy = -np.tanh(1/self.temp)
            exact_m = np.array([0 for i in range(len(self.temp))]) # 정확한 magnetization의 해
        elif self.dim == 2:
            mfa_inst = Newton_Lapson(self.temp, 4)
            mfa_magnets = mfa_inst.get_mfa_magnets()
            mfa_energies = mfa_inst.get_mfa_energies()
            exact_energy = exact_energy_2d(self.temp)
            exact_m = np.power(1-(np.power(1 - np.square(np.tanh(1/self.temp)), 4) / (16 * np.power(np.tanh(1/self.temp), 4))), 1/8)
            exact_m = np.nan_to_num(exact_m)

        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax1.set_title('mean energy')
        #ax1.set_ylim(-2, 0)
        ax1.set_xlabel('T [K]')
        ax1.set_ylabel('Energy per spin')
        for i, data in enumerate(self.data):
            if plot_mode == 'plot':
                ax1.plot(self.temp, data[self.field['energy_m']], self.linestyles[i], label='N =' + data[self.field['size']])
            else:
                ax1.errorbar(self.temp, data[self.field['energy_m']], yerr=data[self.field['std_err_eng']], capsize=3, label='N =' + data[self.field['size']], fmt='.--')
            ax1.legend()
        ax1.plot(self.temp, mfa_energies, label="MFA")
        ax1.legend()
        ax1.plot(self.temp, exact_energy, label="exact solution")
        ax1.legend()

        ax2 = fig.add_subplot(122)
        ax2.set_title('mean magnetization')
        ax2.set_ylim(-0.01, 1)
        ax2.set_xlabel('T [K]')
        ax2.set_ylabel('Magnetization per spin')
        for i, data in enumerate(self.data):
            if plot_mode == 'plot':
                ax2.plot(self.temp, data[self.field['magnet_m']], self.linestyles[i], label='N =' + data[self.field['size']])
            else:
                ax2.errorbar(self.temp, data[self.field['magnet_m']], yerr=data[self.field['std_err_mag']], capsize=3, label='N =' + data[self.field['size']], fmt='.--')
            ax2.legend()
        ax2.plot(self.temp, mfa_magnets, label='MFA')
        ax2.legend()
        ax2.plot(self.temp, exact_m, '-', label="exact solution")
        ax2.legend(loc=(0.72, 0.42))      
        plt.show()

    def show_graph_cx(self, plot_mode):
        '''
        온도에 따른 열용량 및 자화율 그래프 그리는 메소드
        '''
        if self.dim == 1:
            exact_capacity = 1 / ((np.power(np.cosh(1/self.temp), 2)) * (self.temp ** 2))
        elif self.dim == 2:
            q = 2 * np.sinh(2/self.temp) / np.power(np.cosh(2/self.temp), 2)
            phi = np.linspace(0, np.pi/2, len(self.temp))
            K1 = ellipk(q**2)
            E1 = ellipe(q**2)
            exact_capacity = (4/np.pi) * np.power(1/(self.temp*np.tanh(2/self.temp)), 2) * (K1-E1-(1-np.power(np.tanh(2/self.temp), 2))*(np.pi/2 + (2*np.power(np.tanh(2/self.temp), 2)-1)*K1))

        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax1.set_title('Heat capacity')
        ax1.set_xlabel('T [K]')
        ax1.set_ylabel('C/N')
        if self.dim == 1:
            #ax1.set_ylim(0, 1)
            pass
        elif self.dim == 2:
            #ax1.set_xlim(0.5, 4.1) # 2차원의 경우
            #ax1.set_ylim(0, 3) # 2차원의 경우
            ax1.set_ylim(0, 4)
            pass
        for i, data in enumerate(self.data):
            if plot_mode == 'plot':
                ax1.plot(self.temp, data[self.field['capacity']], self.linestyles[i], label='N =' + data[self.field['size']])
            else:
                ax1.errorbar(self.temp, data[self.field['capacity']], yerr=data[self.field['std_err_cap']], capsize=3, fmt='.--', label='N =' + data[self.field['size']])
            ax1.legend()
        ax1.plot(self.temp, exact_capacity, label='exact solution')
        ax1.legend()

        ax2 = fig.add_subplot(122)
        ax2.set_title('magnetic susceptibility')
        ax2.set_xlabel('T [K]')
        ax2.set_ylabel('X/N')
        if self.dim == 1:
            #ax2.set_ylim(0, 18)
            pass
        elif self.dim == 2:
            #ax2.set_ylim(0, 3) # 2차원의 경우
            pass
        for i, data in enumerate(self.data):
            if plot_mode == 'plot':
                ax2.plot(self.temp, data[self.field['suscep']], self.linestyles[i], label='N =' + data[self.field['size']])
            else:
                ax2.errorbar(self.temp, data[self.field['suscep']], yerr=data[self.field['std_err_sus']], capsize=3, fmt='.--', label='N =' + data[self.field['size']])
            ax2.legend()
        plt.show()


#data_inst = ProcessData('extract_quantity.txt', 5)  # 2차원
data_inst = ProcessData('extract_one_d_quantity.txt', 3) # 1차원
data = data_inst.get_data()
field = data_inst.get_field()

#======two_d ising model
#graph_inst = DrawGraph(data, field, 2)
#graph_inst.show_graph('plot')
#graph_inst.show_graph_cx('plot')
#graph_inst.show_graph('errorbar')
#graph_inst.show_graph_cx('errorbar')
#re_data_inst = RearrangeData(data, field, 5)
#re_data_inst.save_re_data('two_d_rearranged_data.txt')
#===========one_d ising model
graph_1d_inst = DrawGraph(data, field, 1)
graph_1d_inst.show_graph('plot')
#graph_1d_inst.show_graph_cx('plot')
graph_1d_inst.show_graph('errorbar')
#graph_1d_inst.show_graph_cx('errorbar')
re_data_1d_inst = RearrangeData(data, field, 3)
re_data_1d_inst.save_re_data('one_d_rearranged_data.txt')