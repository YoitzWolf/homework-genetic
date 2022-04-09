# importing libraries

from typing import Callable, List, Tuple
import numpy as np
import math
from random import choice, randint

from .strAllele import Allele

from .errors import *

class Individual:
    __ALLELES: int # count of Gene Alleles 
    __DIM: int # Size of each gen
    __chromosome: 'np.ndarray[Allele]' # Chromosome.
        #Each gene is integer value with fixed __DIM bit length.
    
    def __str__(self) -> str:
        return f"<Individual: {self.__ALLELES} {self.__DIM}> {list([str(i) for i in self.__chromosome])} ;" 

    def __repr__(self) -> str:
        return self.__str__()

    def get_vector(self) -> np.array:
        w = [
            gene.get_data()
            for gene in self.__chromosome
        ]
        return np.array(w)
    
    def get_allele(self) -> 'np.ndarray[Allele]':
        return self.__chromosome

    def __len__(self) -> int:
        return self.__ALLELES

    def mutate(self, count:int=0, percepts=[1, 1, 1, 1]):
        for i in range(self.__ALLELES):
            for _ in range(count):
                self.__chromosome[i].invert(randint(0, self.__DIM-1))

    def cross(self, other: 'Individual', mutate:float=None) -> 'Individual':
        if len(self) != len(other): raise NOT_EQUAL

        if self.__DIM != len(other.get_allele()[0]) : raise NOT_EQUAL

        child = []

        partner = other.get_allele()

        for i in range(self.__ALLELES):
            child.append(
                Allele(
                    self.__DIM,
                    data=''.join( # !TODO UPDATE THIS PLS
                        [
                            choice(
                                [
                                    self.__chromosome[i].get_bit(j),
                                    partner[i].get_bit(j)
                                ]
                            ) for j in range(self.__DIM)
                        ]
                    )
                )
            )

        #
        # A = [1, 2, 3, 4, 5, 6]
        # B = " ".join(A)
        # B.split(", ")

        child = np.array(child)
        child = Individual(self.__ALLELES, self.__DIM, child)
        if mutate is not None: child.mutate(count=int(self.__DIM * mutate) )

        return child

    def __init__(self, __ALLELES: int, __DIM: int, __chromosome=None, mutate=None):
        if __ALLELES < 2: raise SIZE_ERROR(__ALLELES, "ALLELES")
        if __DIM < 2  : raise SIZE_ERROR(__DIM, "DIM")
        if __chromosome is not None and (len(__chromosome) != __ALLELES or \
            len(__chromosome[0]) != __DIM): raise NOT_EQUAL(__ALLELES)

        self.__ALLELES = __ALLELES
        self.__DIM = __DIM
        # print("DIM ", self.__DIM)
        if __chromosome is None:
            __chromosome = list()
            for i in range(self.__ALLELES):
                __chromosome.append(Allele(self.__DIM))
            __chromosome = np.array(__chromosome)

        self.__chromosome = __chromosome

        if mutate is not None: self.mutate(count=int(self.__DIM * mutate))


class Population():
    __mod: float
    __ACCURACY: int # Accuracy = Allel size
    __DIM: int # count of dimensions of argument (Count of Allels in Individual)
    __MAXES: Tuple[Tuple[float]] # table of maximum values for every dimension
    __population: List[Individual] # List of Individuals
    __function: Callable
    __gen_res: List[Tuple[float]]
    __best: Tuple[Tuple, float]

    def get_accuracy(self):
        return self.__ACCURACY
    
    def get_best(self) -> Tuple[Tuple, float]:
        return self.__best

    def __init__(self, populus_size: int, dim: int, dels: int, maxes: Tuple[float], function: Callable):
        
        self.__ACCURACY = math.ceil(math.log2(dels))# + 1
        
        self.__best = None
        # self.__mod = (max(maxes) - min(maxes)) / dels
        self.__DIM = dim
        self.__MAXES = maxes
        self.__function = function
        self.__gen_res = None
        # print('MAXES', self.__MAXES)
        self.step = 1.0 * ( max(self.__MAXES) - min(self.__MAXES) ) / ((2**(self.__ACCURACY) - 1) + 1)

        self.__population = []
        for _ in range(populus_size):
            self.__population.append(
                Individual(
                    self.__DIM,
                    self.__ACCURACY#,
                    #mutate=0.6
                )
            )
    
    def __modulate(self, x: int) -> 'np.ndarray[Allele]':
        # print("STEP >", f"[ {x} ] ", 2**self.__ACCURACY - 1, x*step)
        return (x * self.step)
        
    def count(self) -> list:
        if self.__gen_res is not None: self.__best = (tuple(), self.__gen_res[0][0])

        self.__gen_res = []
        for i in range(len(self.__population)):
            #print(self.__population[i].get_vector())
            xs = np.array([ self.__modulate(x) for x in self.__population[i].get_vector() ])
            #print(xs)
            z = self.__function(xs)
            self.__gen_res. append(
                (*xs, z) # (1, 2, 3, 4, 12412)
            )
            if self.__best is None or self.__best[1] >= z:
                self.__best = (tuple(xs), z)
        return [] + self.__gen_res
    
    def next(self):
        s = sorted(
            enumerate(self.__population), key=lambda x: self.__gen_res[x[0]][-1]
        )
        # 
        # [a, b, c] -> enumerate([..]) = [(0, a) (1, b) (2, c)]
        #
        # a = lambda x: x
        # (a(y) == y) == True 
        #

        new_population = []
        # Min Max Step -> Min, Min + Step, Min+2Dtep ... MAX
        for i in range(len(self.__population)//2):
            mutmax = int((0.1 + (self.__gen_res[i][-1] / 100)) * 100) # TODO

            # new_population.append(self.__population[s[i][0]])
            new_population.append(s[i][1])
            new_population.append(
                #self.__population[s[i][0]].cross(self.__population[s[i+1][0]], mutate=(randint(0, mutmax)/100) )
                s[i][1].cross(s[i+1][1], mutate=(randint(0, mutmax)/100) )
            )
        
        if len(new_population) > len(self.__population):
            new_population = new_population[:len(self.__population)]
        
        elif len(new_population) < len(self.__population):
            a = choice(self.__population)
            new_population.append(a)

        a = choice(s[len(self.__population)//2:])[1]
        new_population[-1] = a

        self.__population = new_population


'''        
pop = Population(10, 2, 100, ((-10, 10), (-10, 10)), lambda x, y: x+y )
pop.count()
pop.next()
pop.count()
pop.next()
r = pop.count()
print(r)'''