import math
import json
import numpy as np


class BaseModel:
    def __init__(self, conf):
        self.conf = conf

    def utility_function(self,
                         cap,
                         capprima,
                         workS,
                         stocasticproc):
        aux = (stocasticproc*pow(cap, self.conf["alfa"])*pow(workS, 1-self.conf["alfa"]) +
               (1-self.conf[delta])*cap)/(1+self.conf["gamma"])*(1 + self.conf["n"])
        if aux < 0:
            return -10000
        else:
            return math.log(stocasticproc*pow(cap, self.conf["alfa"])*pow(workS, 1-self.conf["alfa"]) -
                          (1+self.conf["gamma"])*(1 + self.conf["n"])*capprima +
                          (1-self.conf[delta])*cap)+self.conf[delta]*math.log(1-workS)

    def calculate_matrix_M(self):
        kvalues = linspace(self.conf["Kmin"], self.conf["Kmax"], self.conf["Kresolution"])
        hvalues = linspace(self.conf["Hmin"], self.conf["Hmax"], self.conf["Hresolution"])
        bvalues = linspace(self.conf["Bmin"], self.conf["Bmax"], self.conf["Bresolution"])
        toret = [[0 for x in range(self.conf["Kresolution"]*self.conf["Bresolution"])] for y in range(0, self.conf["Kresolution"]*self.conf["Hresolution"])]

        countc = 0
        auxc = 0
        for col in range(0, (self.conf["Kresolution"]-1)*(self.conf["Hresolution"]-1)):
            if countc == self.conf["Kresolution"]-1:
                countc = 0
                auxc = auxc + 1
            else:
                countc = countc + 1
            countf = 0
            auxf = 0
            for fila in range(0, (self.conf["Kresolution"]-1)*(self.conf["Bresolution"]-1)):
                if countf == self.conf["Kresolution"]-1:
                    print("reset")
                    countf = 0
                    auxf = auxf + 1
                else:
                    print("++")
                    countf = countf + 1
              
               
                toret[col][fila] = W(kvalues[countc], kvalues[countf], hvalues[auxf], bvalues[auxc])

        return toret


if __name__ == '__main__':
    pass

