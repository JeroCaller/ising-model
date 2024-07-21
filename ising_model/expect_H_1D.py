import matplotlib.pyplot as plt
import numpy as np

class CalcH():
    def __init__(self):
        self.T = np.linspace(0, 10, 100)
        self.N = np.linspace(0, 100, 100)
        self.ex_H = 0

    def calc_T(self):
        self.ex_H = - np.tanh(1/self.T)
        plt.plot(self.T, self.ex_H)
        plt.xlabel('T')
        plt.ylabel('<H>/N', rotation=0)
        plt.show()

    def calc_N(self):
        T = 2
        self.ex_H = - self.N * np.tanh(1/T)
        plt.plot(self.N, self.ex_H)
        plt.show()

    
graph = CalcH()
graph.calc_T()