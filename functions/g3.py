from cmath import pi
import numpy as np
import math as m

from ..genetic.genetic import Individual

def G3(x: np.array, N:int=2):
    return np.sqrt(N)**N * np.prod(x)


def G3Test(x: np.array, N: int=2):
    return ( 1 == sum([y**2 for y in x]) )

def G3Modulator(x: Individual):
    ch = x.get_allele()
    return (ch[0].get_data()*np.cos(ch[1].get_data), ch[0].get_data()*np.sin(ch[1].get_data))

# n = 2
# print(G3(np.array([1/n**0.5, 1/n**0.5])))