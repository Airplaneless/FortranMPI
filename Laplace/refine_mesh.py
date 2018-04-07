import numpy as np
import matplotlib.pylab as plt

if __name__ == '__main__':

    phi = list()
    with open('result', mode='r') as f:
        nx = int(f.readline().strip())
        ny = int(f.readline().strip())
        for line in f.readlines():
            phi_str = line.strip().split(' ')
            phi_ls = list()
            for s in phi_str:
                try:
                    phi_ls.append(float(s))
                except:
                    pass
            phi.append(phi_ls)
    phi = np.array(phi)
    plt.imshow(phi, cmap='hsv')
    plt.colorbar()
    plt.show()