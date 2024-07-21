# 시분초로 나타난 시간을 모두 초로 바꿈

import numpy as np
import matplotlib.pyplot as plt
'''
str_times = ['0:02:22', '0:37:26', '2:44:06', '4:50:42']
time_datas = [data.split(':') for data in str_times]
for i, data in enumerate(time_datas):
    for j in range(3):
        data[j] = int(data[j])

total_sec_list = [data[0] * 3600 + data[1] * 60 + data[2] for data in time_datas]
print(total_sec_list)

exp_sec = [np.format_float_scientific(data, precision=3) for data in total_sec_list]
print(exp_sec)

# ==== chart start ====
lattice_size = [100, 400, 900, 1225]
smallest_digit = len(str(lattice_size[0])) - 1
lattice_size_scale = [data / (10 ** smallest_digit) for data in lattice_size]
print(lattice_size_scale)
smallest_digit_time = len(str(total_sec_list[0])) - 1
total_sec_scale = [data / (10 ** smallest_digit_time) for data in total_sec_list]
print(total_sec_scale)
plt.plot(lattice_size_scale, total_sec_scale, 'o-')
plt.xlabel('lattice size N ' + r'$(10^{})$'.format(smallest_digit))
plt.ylabel('time ' + r'$(10^{} sec)$'.format(smallest_digit_time))
plt.title("Simulation execution time cost")
plt.show()
'''
# ==========
def onlysec(str_times, want_exp=False):
    time_datas = [data.split(':') for data in str_times]
    for i, data in enumerate(time_datas):
        for j in range(3):
            data[j] = int(data[j])
    total_sec_list = [data[0] * 3600 + data[1] * 60 + data[2] for data in time_datas]
    exp_sec = [np.format_float_scientific(data, precision=3) for data in total_sec_list]
    if want_exp is True:
        return total_sec_list, exp_sec
    else:
        return total_sec_list


def showchart(times, lattice_size):
    smallest_digit = len(str(lattice_size[0])) - 1
    lattice_size_scale = [data / (10 ** smallest_digit) for data in lattice_size]
    smallest_digit_time = len(str(times[0])) - 1
    total_sec_scale = [data / (10 ** smallest_digit_time) for data in times]

    plt.plot(lattice_size_scale, total_sec_scale, 'o-', label='data')
    # 대조군
    #plt.plot(lattice_size_scale, np.exp(lattice_size_scale), '^-')
    #const = 1/9.4
    const = 1/275
    temp_var = np.linspace(1, 1000, 1000)
    plt.plot(temp_var, const * np.power(temp_var, 2), '^-', label=r'$t = {:.4}N^2$'.format(const))
    #plt.plot(lattice_size_scale, 9.7 * np.array(lattice_size_scale), '^-')
    # 대조군 끝

    plt.ylim([-10, max(total_sec_scale)+100])
    #plt.xlabel('lattice size N ' + r'$(10^{})$'.format(smallest_digit))
    #plt.ylabel('time t ' + r'$(10^{} sec)$'.format(smallest_digit_time))
    
    plt.xlabel('lattice size N ')
    plt.ylabel('time t ' + '(sec)')

    plt.title("Simulation execution time cost")
    plt.legend()
    plt.show()


times = ['0:00:00', '0:00:01', '0:00:10', '0:00:46', '0:16:04', '1:00:48']
latt_size = [2, 10, 50, 100, 500, 1000]

sec_times = onlysec(times)

print("lattice_size : time_cost")
for i in range(len(sec_times)):
    print("{} : {}".format(latt_size[i], sec_times[i]))

showchart(sec_times, latt_size)