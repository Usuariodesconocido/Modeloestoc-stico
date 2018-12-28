import math
import json
import numpy as np
#Definir la economía
kmin = 0.1
kmax = 4
Hmin = 0.06
Hmax = 0.7
Bmin =1
Bmax =3

delta = 2.5213
alfa = 0.3
gamma = 0.0139
beta = 0.9575
n = 0.0173
sigma = 0.0262

def W(uno, dos, tres, cuatro):
    aux = (cuatro*pow(uno, alfa)*pow(tres, 1-alfa) + (1-delta)*uno)/(1+gamma)*(1 + n)
    if aux < dos:
        return -10000
    else:
        return np.log(cuatro*pow(uno, alfa)*pow(tres, 1-alfa) - (1+gamma)*(1 + n)*dos + (1-delta)*uno)+delta*np.log(1-tres)


lado = 4
colit = 4
fiit = 4

from numpy import linspace
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

for elem in toret:
    print(json.dumps(elem))


V=np.zeros((lado*fiit,1))