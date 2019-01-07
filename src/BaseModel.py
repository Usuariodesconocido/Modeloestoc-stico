import math
import json
import numpy as np
from decimal import Decimal
from scipy import special
import logging
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm


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

        aux = (stocasticproc*cap**float(conf["alfa"])*workS**(1-float(conf["alfa"])) +
               (1-float(conf["sigma"]))*cap)-((1+float(conf["gamma"]))*(1 + float(conf["n"]))*capprima)

        if aux < 0:
            return -10000
        else:

            values = math.log((stocasticproc*(cap**float(conf["alfa"]))*(workS**(1-float(conf["alfa"])))) + ((1-float(conf["sigma"]))*cap) - (
                (1+float(conf["gamma"]))*(1 + float(conf["n"]))*capprima))+(float(conf["delta"])*math.log(float(1-workS)))

            return values

    def calculate_bvalues(self, conf):

        th = (pow(float(conf["ke"]), -1*float(conf["alfa"]))
              * (pow(float(conf["he"]), float(conf["alfa"])-1)))

        mu = (1-float(conf["rho"]))*th
        sig = 2*mu

        s = np.zeros((int(conf["bresolution"]), 1))

        pi = np.zeros((int(conf["bresolution"]), int(conf["bresolution"])))

        s[0, :] = mu/(1-float(conf["rho"]))-int(conf["m"]) * \
            math.sqrt(pow(sig, 2) /

                      (1-(pow(float(conf["rho"]), 2))))
        s[int(conf["bresolution"])-1, :] = mu/(1-float(conf["rho"])) + \
            int(conf["m"])*math.sqrt(pow(sig, 2) /
                                     (1-(pow(float(conf["rho"]), 2))))

        step = (s[int(conf["bresolution"])-1, :] - s[0, :]) / \
            (int(conf["bresolution"])-1)

        for i in range(1, int(conf["bresolution"])-1):
            s[i, :] = s[i-1, :]+step

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

        kvalues = np.linspace(
            float(cc["kmin"]), float(cc["kmax"]), float(cc["kresolution"]))
        hvalues = np.linspace(
            float(cc["hmin"]), float(cc["hmax"]), float(cc["hresolution"]))
        
        bvalues = self.calculate_bvalues(ca)

        toret = [[0 for x in range(int(float(cc["kresolution"])*int(float(cc["bresolution"]))))]
                 for y in range(0, int(float(cc["kresolution"])*int(float(cc["hresolution"]))))]

        countc = -1
        auxc = 0
        for col in range(0, int(float(cc["kresolution"]))*int(float(cc["bresolution"]))):
            if countc == float(cc["kresolution"])-1:
                countc = 0
                auxc = auxc + 1
            else:
                countc = countc + 1

            countf = -1
            auxf = 0
            for fila in range(0, int(float(cc["kresolution"]))*int(float(cc["hresolution"]))):
                if countf == float(cc["kresolution"])-1:
                    countf = 0
                    auxf = auxf + 1
                else:
                    countf = countf + 1

                b = self.utility_function(cc,
                                          kvalues[countc],
                                          kvalues[countf],
                                          hvalues[auxf],
                                          bvalues[auxc])

                toret[col][fila] = b
        return toret


# aldskjfhañdkfhsñdfdsalkñfa´lfdak´lfjal´ksfja´kldfa´lsdfhah


    def autorregresive_matrix(self, conf):

        th = (pow(float(conf["ke"]), -1*float(conf["alfa"]))
              * (pow(float(conf["he"]), float(conf["alfa"])-1)))

        mu = (1-float(conf["rho"]))*th
        sig = 2*mu

        s = np.zeros((int(conf["bresolution"]), 1))

        pi = np.zeros((int(conf["bresolution"]), int(conf["bresolution"])))

        s[0, :] = mu/(1-float(conf["rho"]))-int(conf["m"]) * \
            math.sqrt(pow(sig, 2) /

                      (1-(pow(float(conf["rho"]), 2))))
        s[int(conf["bresolution"])-1, :] = mu/(1-float(conf["rho"])) + \
            int(conf["m"])*math.sqrt(pow(sig, 2) /
                                     (1-(pow(float(conf["rho"]), 2))))

        step = (s[int(conf["bresolution"])-1, :] - s[0, :]) / \
            (int(conf["bresolution"])-1)

        for i in range(1, int(conf["bresolution"])-1):
            s[i, :] = s[i-1, :]+step

        for j in range(0, int(conf["bresolution"])):
            for k in range(0, int(conf["bresolution"])):
                if k == 0:
                    pi[j][k] = 0.5*special.erfc(-1*((s[0, :]-mu-float(
                        conf["rho"])*s[j, :]+step/2)/sig)/math.sqrt(2))

                if k == (int(conf["bresolution"])):
                    pi[j][k] = 1-(0.5*special.erfc(-1*((s[int(conf["bresolution"]), :] -
                                                        mu-float(conf["rho"])*s[j, :]-step/2)/sig)/math.sqrt(2)))
                else:
                    pi[j][k] = 0.5*special.erfc(-1*((s[k, :]-mu-float(conf["rho"])*s[j, :]+step/2)/sig)/math.sqrt(
                        2))-0.5*special.erfc(-1*((s[k, :]-mu-float(conf["rho"])*s[j, :]-step/2)/sig)/math.sqrt(2))
        
        return pi

    def value_function(self, Mconf=None, ARconf=None):
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

        Initialvalue = np.zeros(
            ((int(cc["bresolution"])*int(cc["kresolution"])), ))
        v = Initialvalue
        V = np.ones((int(cc["bresolution"])*int(cc["kresolution"]), ))

        pi = self.autorregresive_matrix(ca)
        s = self.calculate_bvalues(ca)
        O = self.calculate_matrix_M(cc)
        M = np.array(O)
        

        

        W = np.zeros((int(cc["bresolution"]), int(cc["kresolution"])))

        while np.linalg.norm(V-v,2) > float(0.00000000000000000001):
            v = V

            for fila in range(1, int(cc["kresolution"])+1):
                for col in range(1, int(cc["bresolution"])+1):
                    W[fila-1][col-1] = 0
                    
                    for j in range(1, int(cc["kresolution"])+1):

                            
                        
                            W[fila-1][col-1] = W[fila-1, col-1] + \
                                (pi[fila-1, j-1] * \
                                V[(((j-1)*int(cc["kresolution"]))+col)-1, ])

                    

            BELL = np.zeros((int(float(cc["kresolution"]))*int(float(cc["hresolution"])), int(
                float(cc["kresolution"]))*int(float(cc["bresolution"]))))
            countc = 0
            auxc = 1
            for col in range(1, int(float(cc["kresolution"]))*int(float(cc["bresolution"]))+1):
                if countc == float(cc["kresolution"]):
                    countc = 1
                    auxc = auxc + 1
                else:
                    countc = countc + 1
                    

                countf = 0
                auxf = 1
                for fila in range(1, int(float(cc["kresolution"]))*int(float(cc["hresolution"]))+1):
                    if countf == float(cc["kresolution"]):
                        countf = 1
                        auxf = auxf + 1
                    else:
                        countf = countf + 1

                    i1 = (int(cc["kresolution"])*(countc-1)+auxc)-1
                    i2 = (int(cc["kresolution"])*(auxf-1)+countf)-1

                    BELL[i1][i2] = \
                        M[i1, i2]+((float(cc["beta"])*W[countc-1][countf-1]))
                    
                   
                    
            G = np.argmax(BELL, axis=1)
            g=np.array(G)
            V = np.amax(BELL, axis=1)
        
        
        VV = np.zeros((int(float(cc["bresolution"])), int(
            float(cc["kresolution"]))))



        count=0
        for i in range(1, int(float(cc["bresolution"]))+1):
            for j in range(1, int(float(cc["kresolution"]))+1):

                VV[i-1, j-1] = v[((int(float(cc["kresolution"]))*(i-1)+j))-1, ]

                

        # Representar graficamente VV
        

        kvalues = np.linspace(
            float(cc["kmin"]), float(cc["kmax"]), float(cc["kresolution"]))
        hvalues = np.linspace(
            float(cc["hmin"]), float(cc["hmax"]), float(cc["hresolution"]))

        # Funcion de politica para k
        g=np.array((g)) 
        g=g+1
         
        K1 = np.zeros((int(cc["bresolution"]), int(cc["kresolution"])))
        H1 = np.zeros((int(cc["bresolution"]), int(cc["kresolution"])))     
        for i in range(1, int(cc["bresolution"])+1):
            for j in range(1, int(cc["kresolution"])+1):
                imp= int(((i-1)*int(cc["kresolution"]))+j)  

            
                h = math.ceil(g[int(imp)-1,]/(int(cc["kresolution"])))
                s = g[int(imp)-1,]-(int(cc["kresolution"])*(h-1))
                K1[i-1, j-1] = kvalues[s-1]
                H1[i-1, j-1] = hvalues[h-1]
                #print(K1)
                #print("/////////////////////")
                #print(H1)
            

        
        X, Y = np.meshgrid(kvalues, s)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surface = ax.plot_surface(
            X, Y, VV, rstride=1, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False,)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surface2 = ax.plot_surface(
            X, Y, K1, rstride=1, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False,)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        surface3 = ax.plot_surface(
            X, Y, H1, rstride=1, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False,)
        plt.show()

        
if __name__ == '__main__':
    pass
