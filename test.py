from Bplate import *
from sage2d import *

p = 1
c = csi_matrix()
csi = c.final_matrix('file1.dat', 'file1.dat')
s = sage2d()
dod_s1=[]
doa_s1=[]

for i in range(0,1000):
    [beta, f1, f2, CostFunction1, CostFunction2] = s.sage(csi[:, :, 0,i], p)
    print(arcsin(abs(f1) / 0.5)*180/pi)


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