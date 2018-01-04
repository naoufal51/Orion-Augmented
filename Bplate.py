import struct
import numpy as np
from ctypes import *

"""final_matrix

Corrected csi matrix for one sub-carrier.

"""

#
# def final_matrix(csi_corr, csi_target):
#
#     csi_cr = csi_matrix(csi_corr)
#     csi_tr = csi_matrix(csi_target)
#
#     rl_dod = np.zeros(3)
#     rl_doa = np.zeros(3)
#
#     for i in range(0, 3):
#         rl_dod[i] = np.angle(csi_cr[0, i, 0, 0] * np.conjugate(csi_cr[1, i, 0, 0]))
#     for i in range(1, 3):
#         rl_doa[i] = np.angle(csi_cr[0, 0, 0, 0] * np.conjugate(csi_cr[0, i, 0, 0]))
#     csi_tro = np.zeros((2, 3, 1, np.size(csi_tr, 3)), dtype=complex)
#     for k in range(0, np.size(csi_tr, 3) - 1):
#
#         csi_tro[0, :, 0, k] = csi_tr[0, :, 0, k]
#
#         for i in range(0, 3):
#             csi_tro[1, i, 0, k] = csi_tr[1, i, 0, k] * np.exp(1j * rl_dod[i])
#
#         for i in range(1, 3):
#             csi_tro[:, i, 0, k] = csi_tro[:, i, 0, k] * np.exp(1j * rl_doa[i])
#     return csi_tro



"""csi_matrix, reads the csi matrix , removes the spatial mapping,
and stores it in a 'processing friendly' array.

    Args:
        csi: filename of raw csi data file.

    Returns:
        array containing the processed csi matrix.

   """


class csi_matrix:
    def __init__(self):
        pass


    def csi(self,csi):

        triangle = [1, 3, 6]
        broken_perm = 0
        filename = csi
        c = self.read_from_file(filename)
        de_csi = np.transpose(c[0].csi)
        [m, n, s] = de_csi.shape
        csio = np.zeros((s,n,m), dtype=complex)
        num = np.size(c)
        op = np.zeros((m, n, s, num), dtype=np.complex_)
        for i in range(0, np.size(c) - 1):

            perm = c[i].perm
            perm= np.array(perm)
            Nrx = c[i].Nrx
            if Nrx == 1:
                continue
            if np.sum(perm) != triangle[Nrx-1]:
                broken_perm = 1
                print('Invalid CSI File %s with Nrx=%d ',filename,Nrx)
            else:
                c[i].csi=np.array(c[i].csi)
                csio[:, :, :]=(c[i].csi)[:,perm-1,:]
                ret = self.remove_sm(np.transpose(csio), c[i].rate)
                op[:, :, :, i] = ret
        return op

    def remove_sm(self,csi, rate):

        sm_1 = 1

        sm_2_20 = np.matrix([[1, 1], [1, -1]]) / np.sqrt(2)

        sm_2_40 = np.matrix([[1, 1],[1j, 1]]) / np.sqrt(2)

        sm_3_20 =np.matrix( [[-2 * np.pi / 16, -2 * np.pi / (80 / 33), 2 * np.pi / (80 / 3)],
                   [2 * np.pi / (80 / 23), 2 * np.pi / (48 / 13), 2 * np.pi / (240 / 13)],
                   [-2 * np.pi / (80 / 13), 2 * np.pi / (240 / 37), 2 * np.pi / (48 / 13)]])

        sm_3_20 = np.exp(1j * sm_3_20) / np.sqrt(3)


        sm_3_40 = np.matrix([[-2*np.pi/16, -2*np.pi/(80/13), 2*np.pi/(80/23)],
                            [2*np.pi/(80/37), 2*np.pi/(48/11), 2*np.pi/(240/107)],
                            [-2*np.pi/(80/7) ,2*np.pi/(240/83), 2*np.pi/(48/11)]])
        sm_3_40 = np.exp(1j * sm_3_40) / np.sqrt(3)

        m = np.size(csi, 0)
        n = np.size(csi, 1)
        s = np.size(csi, 2)
        ret = np.zeros((m, n, s), dtype=np.complex_)

        if m == 1:
            ret = csi
            return ret
        sm = []
        cond = (np.bitwise_and(rate, 2048) == 2048)

        if cond:
            if m == 3:
                sm = sm_3_40
            elif m == 2:
                sm = sm_2_40
        else:
            if m == 3:
                sm = sm_3_20
            elif m == 2:
                sm = sm_2_20
        for i in range(0, s):
            t = np.array(csi)[:, :, i]
            ret[:, :, i] = np.transpose(np.dot(np.transpose(t),sm.getH()))

        return ret

    #         ret[:, :, i] = np.matrix(np.squeeze(t)).T * np.matrix(sm_2_20).T
    def read_from_file(self,csi_log_file):

        with open(csi_log_file, 'rb') as f:
            data = f.read()

        cur = 0
        length = len(data)
        csi = []

        while cur < (length - 3):
            (field_length,) = struct.unpack(">H", data[cur:cur + 2])
            (code,) = struct.unpack("B", data[cur + 2:cur + 3])
            cur += 3

            csi_bytes = data[cur:cur + field_length - 1]
            cur = cur + field_length - 1

            if code != 187:
                print("Unhandled code %d " % code)
                continue
            iwlnstruct = iwlnl_struct(csi_bytes, from_file=True)
            csi.append(iwlnstruct)
        return csi



class iwlnl_struct:
    def __init__(self, raw_data=False, from_file=False):
        if raw_data:
            self.parse(raw_data, from_file)
    def parse(self, raw_data, from_file=False):
        if not from_file:
            self.unpacked = raw_data[38:]  # skip NETLINK header stuff
        else:
            self.unpacked = raw_data[1:]  # skip the first byte

            # exract noise a b c, bfee_count each held on an unsigned char
            # 0 - 4
            tmp = struct.unpack("BBBBB", self.unpacked[:5])
            self.noise_a = tmp[0]  # 0
            self.noise_b = tmp[1]  # 1
            self.noise_c = tmp[2]  # 2
            self.bfee_count = tmp[3] + (tmp[4] << 8)  # 3
            # extract Nrx, Ntx, rssi_a up to antenna_sel
            # 7 - 19
            tmp = struct.unpack("BBBBBbBBBBBB", self.unpacked[7:19])

            self.Nrx = tmp[0]  # 7
            self.Ntx = tmp[1]  # 8
            self.rssi_a = tmp[2]  # 9
            self.rssi_b = tmp[3]  # 10
            self.rssi_c = tmp[4]  # 11
            self.noise = tmp[5]  # 12
            self.agc = tmp[6]  # 13
            self.antenna_sel = tmp[7]  # 14
            self.length = tmp[8] + (tmp[9] << 8)  # 15-16
            self.rate = tmp[10] + (tmp[11] << 8)  # 16-17
            # number of subcarriers
            self.nrs = 30
            self.perm = []
            self.perm.append(((self.antenna_sel) & 0x3) + 1)
            self.perm.append(((self.antenna_sel >> 2) & 0x3) + 1)
            self.perm.append(((self.antenna_sel >> 4) & 0x3) + 1)
            # print self.perm
            self.csi = self.parse_csi(self.unpacked[19:])

        # """ Given raw_data (bytes) parse CSI complex values """

    def parse_csi(self, raw_data):
        index = 0
        remainder = 0
        # make a list of 30 elements
        csi = [None] * 30
        try:
            for i in range(0, self.nrs):
                index += 3
                remainder = (index % 8)
                Hx = np.matrix(np.zeros((self.Nrx, self.Ntx), complex))
                for r in range(0, self.Nrx):
                    for t in range(0, self.Ntx):
                        #            first = struct.unpack("B",raw_data[index/8])[0] >> remainder
                        first = struct.unpack('B', bytes([raw_data[index // 8]]))[0] >> remainder
                        second = (struct.unpack('B', bytes([raw_data[index // 8 + 1]]))[0] << (8 - remainder))
                        tmp = (c_byte(first | second).value)
                        real = (c_double(tmp).value)
                        first = (struct.unpack('B', bytes([raw_data[index // 8 + 1]]))[0] >> remainder)
                        second = (struct.unpack('B', bytes([raw_data[index // 8 + 2]]))[0] << (8 - remainder))
                        tmp = (c_byte(first | second).value)
                        imag = (c_double(tmp).value)
                        index += 16
                        Hx.itemset((r, t), complex(real, imag))
                csi[i] = Hx
        except IndexError:
            pass
        return csi