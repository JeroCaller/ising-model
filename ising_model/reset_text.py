import os


def reset_text(file_name):
    try:
        with open(file_name, 'w') as f:
            pass
    except Exception as e:
        print(e)
    else:
        print('data successfully deleted.')


def read_txt(file_name):
    try:
        with open(file_name, 'r') as f:
            data = f.readlines()
    except Exception as e:
        print(e)
    return data


def delete_file(file_name):
    try:
        os.remove(file_name)
    except Exception as e:
        print(e)
    else:
        print('The file is succefully deleted.')


#delete_file('two_d_quantity.txt')
reset_text('two_d_quantity.txt')
#print(read_txt('aaa.txt'))