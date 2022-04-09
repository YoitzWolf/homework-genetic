from cmath import pi
import numpy as np
import math as m


def Zakharov(x: np.array, N:int=2):
    rez = 0

    rez += sum (
        [
           (x[i]*x[i]) for i in range(N)
        ]
    )

    rez += sum(
        [
            x[i]*(i+1)*0.5
            for i in range(N)
        ]
    )**2

    rez += sum(
        [
            x[i]*(i+1)*0.5
            for i in range(N)
        ]
    )**4


    return rez