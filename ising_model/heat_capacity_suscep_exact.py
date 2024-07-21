'''
1, 2차원의 스핀 당 열용량과 자화율의 정확한 해를 온도에 따라 그래프로 그려본다.
'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipk
from scipy.special import ellipe
from scipy.integrate import quad
import math

def one_capacity(ts):
    # ts == list of temperatures
    C = 1 / ((np.power(np.cosh(1/ts), 2)) * (ts ** 2))
    return C

def one_suscep(ts):
    X = np.array([0 for t in range(len(ts))])
    return X

def two_capacity(ts):
    q = 2 * np.sinh(2/ts) / np.power(np.cosh(2/ts), 2)
    phi = np.linspace(0, np.pi/2, len(ts))
    K1 = ellipk(q**2)
    #E1 = ellipe((q**2)*(np.power(np.sinh(phi), 2))/(np.power(np.sin(phi), 2)))
    E1 = ellipe(q**2)
    '''def calc_e(x, one_q):
        return np.sqrt(1-(one_q**2)*np.power(np.sinh(x), 2))'''

    #E2 = np.array([quad(calc_e, 0, np.pi/2, args=(q[i],))[0] for i in range(len(ts))])
    C = (4/np.pi) * np.power(1/(ts*np.tanh(2/ts)), 2) * (K1-E1-(1-np.power(np.tanh(2/ts), 2))*(np.pi/2 + (2*np.power(np.tanh(2/ts), 2)-1)*K1))
    return C

def two_capacity_app(ts):
    # 책참조
    Tc = 2 / np.log(1+np.sqrt(2))
    C = -(2/np.pi)*np.power(2/Tc, 2) * np.log(np.abs(1-ts/Tc))
    return C

def two_capacity_app2(ts):
    # 블로그 참조
    Tc = 2 / np.log(1+np.sqrt(2))
    C = (2/np.pi)*np.power(2/Tc, 2) * (-np.log(np.abs(1-ts/Tc)) + np.log(Tc/2) - (1 + np.pi/4) )
    return C


def two_suscep(ts):
    '''
    정확한 해라기 보단 근사적으로 푼 해
    '''
    Tc = 2 / np.log(1+np.sqrt(2))
    X = np.power(np.abs(ts-Tc), -7/4)
    #X = (ts-Tc) ** (7/4)
    #X = 2 ** (7/4)
    return X / 400
    #return X
    

def draw_graph(ts, quantity, title, q_name):
    plt.xlabel('T [K]')
    plt.ylabel(q_name)
    plt.title(title)
    #plt.xlim(0, 10)
    #plt.ylim(0, 10)
    plt.plot(ts, quantity)
    plt.show()


if __name__ == '__main__':
    temps = np.arange(0.1, 4.0 + 0.1, 0.1)
    one_c = one_capacity(temps)
    one_x = one_suscep(temps)
    two_c = two_capacity(temps)
    two_c_app = two_capacity_app(temps)
    two_c_app2 = two_capacity_app2(temps)
    two_x = two_suscep(temps)

    titles = ['one dimensional specific heat capacity per spin', 'one dimensional zero-field susceptibility per spin',
    'two dimensional capacity per spin', 'two dimensional zero-field susceptibility per spin']
    data_set = [one_c, one_x, two_c_app, two_x]
    q_names = ['specific heat capacity per spin', 'zero-field susceptibility per spin']

    draw_graph(temps, two_x, titles[3], q_names[0])