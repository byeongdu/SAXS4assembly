import h5py
import numpy as np

def read(filename):
    h5 = h5py.File(filename, 'r')
    dt = h5.get('dt')
    dt = np.array(dt)
    return dt

def saveas(filename, res):
    hf = h5py.File(filename, 'w')
    hf.create_dataset('dt', data=res)
    hf.close()