from cmath import pi
import numpy as np
import math as m

# from ..genetic.genetic import Individual

def G3(x: np.array, N:int=2):
    rs = x[0] * x[1]# / x[0]
    # for i in range(1, N):
    #     rs *= x[i]
    #print(rs, x[0])
    #input()
    return rs * ((np.sqrt(N))**N)

def G3Test(x: np.array, N: int=2):
    return ( 1 == sum([y**2 for y in x]) )

def G3Modulator(x):
    ch = x.get_allele()
    angle = ch[0].get_data()
    ch[1].clear()
    return (np.cos(angle), np.sin(angle))

# n = 2
# print(G3(np.array([1/n**0.5, 1/n**0.5])))