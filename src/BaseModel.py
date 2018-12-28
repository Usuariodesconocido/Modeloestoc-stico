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
                         conf,
                         cap,
                         capprima,
                         workS,
                         stocasticproc):
        aux = (stocasticproc*pow(cap, float(conf["alfa"]))*pow(workS, 1-float(conf["alfa"])) +
               (1-float(conf["delta"]))*cap)/(1+float(conf["gamma"]))*(1 + float(conf["n"]))
        if aux < 0:
            return -10000
        else:
            return math.log(stocasticproc*pow(cap, float(conf["alfa"]))*pow(workS, 1-float(conf["alfa"])) -
                            (1+float(conf["gamma"]))*(1 + float(conf["n"]))*capprima +
                            (1-float(conf["delta"]))*cap)+float(conf["delta"])*math.log(1-workS)

    def calculate_matrix_M(self, Mconf=None):
        if Mconf is None:
            if self.Mconf is not None:
                cc = self.Mconf
            else:
                logging.critical("NO CONF")
                return
        else:
            cc = Mconf

        print json.dumps(cc, indent=2)

        kvalues = np.linspace(
            float(cc["kmin"]), float(cc["kmax"]), float(cc["kresolution"]))
        hvalues = np.linspace(
            float(cc["hmin"]), float(cc["hmax"]), float(cc["hresolution"]))
        bvalues = np.linspace(
            float(cc["bmin"]), float(cc["bmax"]), float(cc["bresolution"]))
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
