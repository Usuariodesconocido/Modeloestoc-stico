import math
import json
import numpy as np


class BaseModel:
    def __init__(self, conf):
        self.conf = conf

    def utility_function(self,
                         cap,
                         capprima,
                         tres,
                         cuatro):
        aux = (cuatro*pow(uno, alfa)*pow(tres, 1-alfa) +
               (1-delta)*uno)/(1+gamma)*(1 + n)
        if aux < 0:
            return -10000
        else:
            return np.log(cuatro*pow(uno, alfa)*pow(tres, 1-alfa) -
                          (1+gamma)*(1 + n)*dos +
                          (1-delta)*uno)+delta*np.log(1-tres)

    def calculate_matrix_M(self):
        kvalues = linspace(kmin, kmax, lado)
        hvalues = linspace(Hmin, Hmax, colit)
        tvalues = linspace(Bmin, Bmax, fiit)
        toret = [[0 for x in range(lado*fiit)] for y in range(0, lado*colit)]

        countc = 0
        auxc = 0
        for col in range(0, (lado-1)*(colit-1)):
            if countc == lado-1:
                countc = 0
                auxc = auxc + 1
            else:
                countc = countc + 1
            countf = 0
            auxf = 0
            for fila in range(0, (lado-1)*(fiit-1)):
                if countf == lado-1:
                    print("reset")
                    countf = 0
                    auxf = auxf + 1
                else:
                    print("++")
                    countf = countf + 1
              
               
                toret[col][fila] = W(kvalues[countc], kvalues[countf], hvalues[auxf], tvalues[auxc])

        return toret


if __name__ == '__main__':
    pass

