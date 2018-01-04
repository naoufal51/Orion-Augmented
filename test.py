from Bplate import *
from sage2d import *

p = 1
c = csi_matrix()
csi = c.final_matrix('file1.dat', 'file1.dat')
s = sage2d()
result = s.sage(csi[:, :, :, 1], p)



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