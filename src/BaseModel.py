import math
import json
import numpy as np


class BaseModel:
    def __init__(self,
                 Mconf=None,
                 ARconf=None):
        self.Mconf = Mconf
        self.ARconf = ARconf

    def utility_function(self,
                         conf
                         cap,
                         capprima,
                         workS,
                         stocasticproc):
        aux = (stocasticproc*pow(cap, conf["alfa"])*pow(workS, 1-conf["alfa"]) +
               (1-conf[delta])*cap)/(1+conf["gamma"])*(1 + conf["n"])
        if aux < 0:
            return -10000
        else:
            return math.log(stocasticproc*pow(cap, conf["alfa"])*pow(workS, 1-conf["alfa"]) -
                            (1+conf["gamma"])*(1 + conf["n"])*capprima +
                            (1-conf[delta])*cap)+conf[delta]*math.log(1-workS)

    def calculate_matrix_M(self, Mconf=None):
        if Mconf is None:
            if self.Mconf is not None:
                cc = self.Mconf
            else:
                logging.critical("NO CONF")
                return
        else:
            cc = Mconf
            
        kvalues = linspace(
            cc["Kmin"], cc["Kmax"], cc["Kresolution"])
        hvalues = linspace(
            cc["Hmin"], cc["Hmax"], cc["Hresolution"])
        bvalues = linspace(
            cc["Bmin"], cc["Bmax"], cc["Bresolution"])
        toret = [[0 for x in range(cc["Kresolution"]*cc["Bresolution"])]
                 for y in range(0, cc["Kresolution"]*cc["Hresolution"])]

        countc = 0
        auxc = 0
        for col in range(0, (cc["Kresolution"]-1)*(cc["Hresolution"]-1)):
            if countc == cc["Kresolution"]-1:
                countc = 0
                auxc = auxc + 1
            else:
                countc = countc + 1
            countf = 0
            auxf = 0
            for fila in range(0, (cc["Kresolution"]-1)*(cc["Bresolution"]-1)):
                if countf == cc["Kresolution"]-1:
                    print("reset")
                    countf = 0
                    auxf = auxf + 1
                else:
                    print("++")
                    countf = countf + 1

                toret[col][fila] = utility_function(cc,
                                                    kvalues[countc],
                                                    kvalues[countf],
                                                    hvalues[auxf],
                                                    bvalues[auxc])

        return toret


if __name__ == '__main__':
    pass
