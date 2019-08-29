import numpy as np
import csv


def check_csv(file):
    data = []
    answer = True
    with open(file, 'r') as data_file:
        data = data_file.read()
        if ',' in data:
            print('csv')
            answer = True
        else:
            print('non csv')
            answer = False
    return answer


def with_csv(file):
    data = []
    with open(file, 'r') as data_file:
        print('cvs-collecting')
        csvReader = csv.reader(data_file)
        data = [[float(line) for line in lines] for lines in csvReader]
        data = np.array(data)
    return data

def with_np(file):
    print('np calc')
    data = np.loadtxt(file)
    return data


def data_in(file):
    data = []
    check = check_csv(file)
    if check:
        data = with_csv(file)
    elif not check:
        data = with_np(file)
    else:
        print('error check should be bool')
    data = np.transpose(data)
    return data

def test1():
    '''
    test to check if data in is correct
    '''
    wanted = [0.0123]# i think
    file = 'BslA_NR_WT_ellipsoids/29543_44.dat'
    so = ''
    try:
        data = data_in(file)
        print(data)
    except Exception as e:
        print('error:', e)
    try:
        check, comp_len = compare(data, wanted)
    except Exception as e:
        print('error:', e)
        check = False
    if check:
        so = 'done'
    else:
        so = 'error: comparision failed'
    print(so)

def compare(a, b):
    '''
    used to test
    '''
    minim = min(len(a), len(b))
    same = True
    for i in range(minim):
        same = (a[i] == b[i]) and same
    return same, minim

def test2():
    '''
    test to check check_csv works
    '''
    file_a = 'BslA_NR_WT_ellipsoids/29543_44.dat'
    file_b = 'dppc_pxg/33890_91.dat'
    #check_csv(file_a)
    #check_csv(file_b)
    result = type(data_in(file_a))==type(data_in(file_b))
    print(result)
#test2()
#out = data_in("d2o/29553_54.dat")
#print(out)