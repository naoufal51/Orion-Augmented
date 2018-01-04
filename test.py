from Bplate import *
from sage2d import *

p = 1
c = csi_matrix()
csi = c.final_matrix('file1.dat', 'file1.dat')
s = sage2d()
[beta, f1, f2, CostFunction1, CostFunction2] = s.sage(csi[:, :, 0,0], p)
print(beta,f2, f1)


    # def main():
    #     p = 1
    #     c = csi_matrix()
    #     csi = c.final_matrix('file1.dat', 'file1.dat')
    #     s = sage2d()
    #     length = csi.shape
    #     for i in range(length[3]):
    #         result = s.sage(csi[:, :, :, i], p)
    #     return result
    #
    #
    # if __name__ == '__main__':
    #     main()