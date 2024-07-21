# mean field approx 해를 온도에 따라 수치적으로 구하는 과정. 

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

class Newton_Lapson():
    def __init__(self, temps, nns):
        self.T = temps
        self.one_temp = 1
        self.nns = nns

    def f(self, m):
        formula = np.tanh((self.nns/self.one_temp)*m) - m
        return formula

    def get_many_solutions(self):
        solutions = []
        for t in self.T:
            self.one_temp = t
            solutions.append(fsolve(self.f, [1, 100])[0])
        return solutions

    def get_temps(self):
        return self.T


class DrawGraph():
    def __init__(self, temps, quantity):
        self.T = temps
        self.quantity = quantity

    def draw_mag(self):
        plt.plot(self.T, self.quantity)
        plt.xlabel('temperature [K]')
        plt.ylabel('magnetization m')
        plt.show()

    def draw_ener(self):
        plt.plot(self.T, self.quantity)
        plt.xlabel('temperature [K]')
        plt.ylabel('energy <H>')
        plt.show()


def mft_energy(temps_list, nns, m_list):
    energy = - (nns/2) * np.power(m, 2)
    return energy


temps = np.arange(0.1, 4.0 + 0.1, 0.1)
nns = 2

method_inst = Newton_Lapson(temps, nns)
m = method_inst.get_many_solutions()

energies = mft_energy(temps, nns, m)

graph_inst = DrawGraph(temps, energies)
graph_inst.draw_mag()