import math
import json
import numpy as np
from decimal import Decimal
from scipy import special
import logging

class BaseModel:
    def __init__(self,
                 Mconf=None,
                 ARconf=None):
        self.Mconf = Mconf
        self.ARconf = ARconf

    def utility_function(self,
                         conf,
                         cap,
                         capprima,
                         workS,
                         stocasticproc):
        aux = (stocasticproc*pow(cap, float(conf["alfa"]))*pow(workS, 1-float(conf["alfa"])) +
               (1-float(conf["delta"]))*cap)/(1+float(conf["gamma"]))*(1 + float(conf["n"]))
        if aux < capprima:
            return -10000
        else:
            return math.log1p(Decimal(stocasticproc*pow(cap, float(conf["alfa"]))*pow(workS, 1-float(conf["alfa"])) -
                                      (1+float(conf["gamma"]))*(1 + float(conf["n"]))*capprima +
                                      (1-float(conf["delta"]))*cap))+(float(conf["delta"])*math.log1p(float(1-workS)))

    def calculate_bvalues(self, conf):
        m=3
        s = np.zeros((int(conf["bresolution"]), 1))
        pi = np.zeros((int(conf["bresolution"]), int(conf["bresolution"])))
        s[0, :] = float(conf["mu"])/(1-float(conf["rho"]))-m * \
            math.sqrt(pow(float(conf["sig"]), 2) /
                      (1-(pow(float(conf["rho"]), 2))))
        s[int(conf["bresolution"])-1, :] = float(conf["mu"])/(1-float(conf["rho"])) + \
            m*math.sqrt(pow(float(conf["sig"]), 2) /
                        (1-(pow(float(conf["rho"]), 2))))
        step = (s[int(conf["bresolution"])-1, :] - s[0, :]) / \
            (int(conf["bresolution"])-1)

        for i in range(1, int(conf["bresolution"])-1):
            s[i, :] = s[i-1, :]+step

        for j in range(0, int(conf["bresolution"])-1):
            for k in range(0, int(conf["bresolution"])-1):
                if k == 0:
                    pi[j][k] = 0.5*special.erfc((-1*(s[0, :]-float(conf["mu"])-float(
                        conf["rho"]*s[j, :])+step/2)/float(conf["sig"]))/math.sqrt(2))
                elif k == int(conf["bresolution"]):
                    pi[j][k] = 1-(0.5*special.erfc((-(s[int(conf["bresolution"])-1, :]-float(
                        conf["mu"])-float(conf["rho"]*s[j, :])+step/2)/float(conf["sig"]))/math.sqrt(2)))
                else:
                    pi[j][k] = 0.5*special.erfc((-(s[k, :]-float(conf["mu"])-float(conf["rho"]*s[j, :])+step/2)/float(conf["sig"]))/math.sqrt(
                        2))-0.5*special.erfc((-(s[k, :]-float(conf["mu"])-float(conf["rho"]*s[j, :])+step/2)/float(conf["sig"]))/math.sqrt(2))

        return s

    def calculate_matrix_M(self, Mconf=None, ARconf=None):
        if Mconf is None:
            if self.Mconf is not None:
                cc = self.Mconf
            else:
                logging.critical("NO CONF")
                return
        else:
            cc = Mconf

        if ARconf is None:
            if self.ARconf is not None:
                ca = self.ARconf
            else:
                logging.critical("NO CONF")
                return
        else:
            ca = ARconf

        print(ca)

        kvalues = np.linspace(
            float(cc["kmin"]), float(cc["kmax"]), float(cc["kresolution"]))
        hvalues = np.linspace(
            float(cc["hmin"]), float(cc["hmax"]), float(cc["hresolution"]))

        bvalues = self.calculate_bvalues(ca)

        toret = [[0 for x in range(int(float(cc["kresolution"])*int(float(cc["bresolution"]))))]
                 for y in range(0, int(float(cc["kresolution"])*int(float(cc["hresolution"]))))]

        countc = 0
        auxc = 0
        for col in range(0, int(float(cc["kresolution"])-1)*int(float(cc["hresolution"])-1)):
            if countc == float(cc["kresolution"])-1:
                countc = 0
                auxc = auxc + 1
            else:
                countc = countc + 1
            countf = 0
            auxf = 0
            for fila in range(0, int(float(cc["kresolution"])-1)*int(float(cc["bresolution"])-1)):
                if countf == float(cc["kresolution"])-1:
                    countf = 0
                    auxf = auxf + 1
                else:
                    countf = countf + 1

                toret[col][fila] = self.utility_function(cc,
                                                         kvalues[countc],
                                                         kvalues[countf],
                                                         hvalues[auxf],
                                                         bvalues[auxc])

        return toret


if __name__ == '__main__':
    pass
