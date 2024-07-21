from asyncio.proactor_events import _ProactorBaseWritePipeTransport
import matplotlib.pyplot as plt
import numpy as np

class ProcessData():
    def __init__(self, file_name):
        self.file_name = file_name
        
        self.raw_data = self.read_data()
        self.field = self.extract_field()
        self.all_data = []
        for i in range(len(self.raw_data)):
            one_data = self.process_one_data(self.raw_data[i])
            self.all_data.append(one_data)

    def read_data(self):
        try:
            with open(self.file_name, 'r') as f:
                raw_data = f.readlines()
        except FileNotFoundError:
            print('에러: "{}라는 이름의 파일이 존재하지 않습니다!"'.format(self.file_name))
        else:
            return raw_data

    def extract_field(self):
        field = self.raw_data.pop(0).strip().split(',')
        return field

    def process_one_data(self, one_data=''):
        split_data = one_data.split('@')

        def erase_char(some_property=''):
            new_property = some_property.replace('[', '').replace(']', '').split(' ')
            for i, data in enumerate(new_property):
                if data == '':
                    del new_property[i]
            for i, data in enumerate(new_property):
                new_property[i] = float(data)
            return new_property

        T = erase_char(split_data[0])
        init_latt = erase_char(split_data[1])
        all_latt = split_data[2].replace('{', '').replace('}', '').replace('\'', '').split(',')
        for i, val in enumerate(all_latt):
            all_latt[i] = val.strip().split(': ')
        for i, val in enumerate(all_latt):
            all_latt[i][0] = float(val[0])
            all_latt[i][1] = val[1].replace('[[', '').replace(']]', '').split(' ')
            for j, val2 in enumerate(all_latt[i][1]):
                if val2 == '':
                    del all_latt[i][1][j]
            for j, val2 in enumerate(all_latt[i][1]):
                all_latt[i][1][j] = float(val2)
        for i, val in enumerate(all_latt):
            all_latt[i][1] = np.array(val[1]).reshape(2, 2)
        return [T, init_latt, all_latt]

    def get_data(self):
        return self.all_data


class ShowData():
    def __init__(self, data):
        self.T = data[0]
        self.init_latt = np.array(data[1]).reshape(2, 2)
        self.all_latt = data[2]

        self.sample_T = [0.1, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
        self.sample_latt = []

        self.set_sample_latt()
        print(self.sample_latt)
        ok = self.sample_latt[0]
        not_ok = self.sample_latt[1]
        print(ok)
        print(not_ok)

    def set_sample_latt(self):
        for i, data in enumerate(self.all_latt):
            if data[0] in self.sample_T:
                self.sample_latt.append(data[1])

    def show_graph(self):
        fig = plt.figure()
        print(fig, type(fig))
        cmap = plt.get_cmap('binary')
        
        ax_init = fig.add_subplot(241)
        ax_init.matshow(self.init_latt, cmap=cmap, vmin=-1, vmax=1)
        ax = []
        for i in range(len(self.sample_T)):
            ax.append(fig.add_subplot(2, 4, i+2))
            ax[i].matshow(self.sample_latt[i], cmap=cmap, vmin=-1, vmax=1)

        #plt.clim(-1, 1)
        plt.show()
        


FILE_NAME = 'two_d_lattice.txt'
data_inst = ProcessData(FILE_NAME)
new_data = data_inst.get_data()
#print(new_data)

graph_inst = ShowData(new_data[0])
graph_inst.show_graph()