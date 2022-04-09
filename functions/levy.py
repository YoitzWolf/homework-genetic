from cmath import pi
import numpy as np
import math as m

def Levy(w: np.array, N:int=2):
    return np.sin(pi*w[0])**2 + sum(
            [
                (w_ - 1)**2 * (1 + 10 * np.sin(np.pi*w_ + 1)**2)
                for w_ in w[:-1]
            ]
        ) + (w[N-1]-1)**2 * (1 + np.sin(2*pi*w[N-1])**2 )


def Levy_arg_w(x: float) -> float:
    
    return 1 + (x - 1)/4


def LevyFull(x: np.array, N:int=2):
    w = Levy_arg_w(x)
    return Levy(w, N)