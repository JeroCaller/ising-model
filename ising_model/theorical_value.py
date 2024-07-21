'''
이론값의 그래프 그리기
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk


class DrawTheory():
    def __init__(self, q, T):
        self.q = q
        self.T = T

        self.m = self.mean_field2()
        self.energy = self.mean_field_energy2()

    def mean_field(self):
        '''
        여기서 q는 특정 스핀에 대해서 가장 가까운 스핀의 개수이다.
        1차원: 2, 2차원: 4, 3차원: 6
        볼츠만 상수 및 J는 1로 둔다.
        '''
        Tc = self.q
        m_plus = np.sqrt(3) * np.power(self.T/Tc, 3/2) * np.power((Tc/self.T)-1, 1/2)
        return m_plus

    def mean_field2(self):
        # 이 식이 정확한 듯 하다.
        Tc = self.q
        # m = np.power(Tc-T, 1/2) / np.sqrt(Tc)
        m = np.power(1-(self.T/Tc), 1/2)
        return m

    def mean_field_energy(self):
        Tc = self.q
        energy_expectation = - self.q * (1- self.T/Tc)
        return energy_expectation

    def mean_field_energy2(self):
        # 이 식이 더 정확한 듯 하다.
        Tc = self.q
        energy_expectation = - Tc * np.power(np.tanh(np.sqrt((Tc/self.T)*(Tc/self.T-1))),2)
        return energy_expectation

    def draw_graph(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax1.set_xlabel('T [K]')
        ax1.set_ylabel('<E>')
        ax1.set_title('MFA <E>')
        #ax1.set_xlim(0, 4.0)
        #ax1.set_ylim(-2, 0)
        ax1.plot(self.T, self.energy)
        
        ax2 = fig.add_subplot(122)
        ax2.set_xlabel('T [K]')
        ax2.set_ylabel('m')
        ax2.set_title('MFA magnetization')
        #ax2.set_xlim(0, 4.0)
        #ax2.set_ylim(0, 1)
        ax2.plot(self.T, self.m)
        
        plt.show()


class DrawEnergy2D():
    def __init__(self, T):
        self.T = T

        self.exact_energy_2d_1 = self.calc_energy_2d()
        self.exact_energy_2d_2 = self.calc_energy_2d_ver2()
        self.exact_energy_2d_3 = self.calc_energy_2d_ver3()
        self.exact_energy_2d_4 = self.calc_energy_2d_ver4()
        self.exact_energy_2d_5 = self.calc_energy_2d_ver5()
        self.exact_energy_2d_6 = self.calc_energy_2d_ver6()

    def calc_energy_2d(self):
        q = 2 * np.sinh(2/self.T) / np.power(np.cosh(2/self.T), 2)
        K1_approx = (np.sqrt(2) / (2 * (1 - q**2))) * (np.arctanh(1/np.sqrt(2)) - q * np.arctanh(q / np.sqrt(2)))
        part1 = 1 / np.tanh(2/self.T)
        part2 = 2 * np.power(np.tanh(2/self.T), 2) - 1
        energy = - part1 * (1 + (2/np.pi) * part2 * K1_approx)
        return energy

    def calc_energy_2d_ver2(self):
        '''
        인터넷 계산 사이트에서 얻어옴.
        '''
        q = 2 * np.sinh(2/self.T) / np.power(np.cosh(2/self.T), 2)
        K1_approx = np.pi / 2 + (q**2)*(np.pi ** 3) / 48 + (q**2) * (9 * (q ** 2) - 4) * (np.pi ** 5) / 3840
        part1 = 1 / np.tanh(2/self.T)
        part2 = 2 * np.power(np.tanh(2/self.T), 2) - 1
        energy = - part1 * (1 + (2/np.pi) * part2 * K1_approx)
        return energy

    def calc_energy_2d_ver3(self):
        q = 2 * np.sinh(2/self.T) / np.power(np.cosh(2/self.T), 2)
        K1_approx = (np.pi / 2) * (1 + (q**2)/4)
        part1 = 1 / np.tanh(2/self.T)
        part2 = 2 * np.power(np.tanh(2/self.T), 2) - 1
        energy = - part1 * (1 + (2/np.pi) * part2 * K1_approx)
        return energy

    def calc_energy_2d_ver4(self):
        q = 2 * np.sinh(2/self.T) / np.power(np.cosh(2/self.T), 2)
        K1_approx = np.pi / 2 + (q**2)*(np.pi**3)/48 - (q**2) * (np.pi**5)/960
        part1 = 1 / np.tanh(2/self.T)
        part2 = 2 * np.power(np.tanh(2/self.T), 2) - 1
        energy = - part1 * (1 + (2/np.pi) * part2 * K1_approx)
        return energy

    def calc_energy_2d_ver5(self):
        q = 2 * np.sinh(2/self.T) / np.power(np.cosh(2/self.T), 2)
        K1_approx = np.pi / 2 + (q**2)*(np.pi**3)/48 - (q**2) * (np.pi**5)/960 + (q**2)*(np.pi**7)/40320 - (q**2)*(np.pi**9)/2903040
        part1 = 1 / np.tanh(2/self.T)
        part2 = 2 * np.power(np.tanh(2/self.T), 2) - 1
        energy = - part1 * (1 + (2/np.pi) * part2 * K1_approx)
        return energy

    def calc_energy_2d_ver6(self):
        q = 2 * np.sinh(2/self.T) / np.power(np.cosh(2/self.T), 2)
        K1_approx = ellipk(q**2)
        part1 = 1 / np.tanh(2/self.T)
        part2 = 2 * np.power(np.tanh(2/self.T), 2) - 1
        energy = - part1 * (1 + (2/np.pi) * part2 * K1_approx)
        #energy = np.where(energy <= -2, -2, energy)
        return energy

    def show_graph(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('T [K]')
        ax.set_ylabel('<E>')
        ax.set_title('exact <E> in 2-D')
        ax.set_xlim(0, 4.0)
        ax.set_ylim(-2, 0)
        #plt.xticks([0, 1, 2, 3, 4])
        #plt.yticks([-2, -1, 0])
        #ax.plot(self.T, self.exact_energy_2d_1, label='ver1')
        #ax.legend()
        #ax.plot(self.T, self.exact_energy_2d_2, label='ver2')
        #ax.legend()
        #ax.plot(self.T, self.exact_energy_2d_3, label='ver3')
        #ax.legend()
        #ax.plot(self.T, self.exact_energy_2d_4, label='ver4')
        #ax.legend()
        #ax.plot(self.T, self.exact_energy_2d_5, label='ver5')
        #ax.legend()
        #ax.plot(self.T, self.exact_energy_2d_6, label='ver6')
        #ax.legend()
        ax.plot(self.T, self.exact_energy_2d_6, label='ver6')

        plt.show()


T = np.arange(0.1, 4.0 + 0.1, 0.1)

inst_1d = DrawTheory(2, T)
inst_1d.draw_graph()

#inst_2d = DrawTheory(4, T)
#inst_2d.draw_graph()

#exact_energy_2d = DrawEnergy2D(T)
#exact_energy_2d.show_graph()


#==========
'''
some_q = np.arange(0.1, 1.0 + 0.1, 0.1)
K1 = np.arcsin(some_q) / some_q
plt.plot(some_q, K1)
plt.show()'''