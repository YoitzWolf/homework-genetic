# importing libraries

from typing import Callable, List, Tuple
import numpy as np
import math
from random import choice, randint

from .strAllele import Allele

from .errors import *

def RANDOMCHROMO(first: 'Individual', second: 'Individual'):
    return list(
        [
            list([
                choice(
                    [
                        first[i].get_bit(j),
                        second[i].get_bit(j)
                    ]
                ) for j in range(first.get_allele_size())
            ])
        for i in range(len(first)) ] 
    )


class Individual:
    __ALLELES: int # count of Gene Alleles 
    __DIM: int # Size of each gen
    __chromosome: 'np.ndarray[Allele]' # Chromosome.
        #Each gene is integer value with fixed __DIM bit length.
    __GETCHILDCHROMO: Callable

    def __getitem__(self, arg):
        return self.__chromosome[arg]

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
    
    def get_list(self) -> list:
        w = [
            gene.get_data()
            for gene in self.__chromosome
        ]
        return list(w)
    
    def get_allele(self) -> 'np.ndarray[Allele]':
        return self.__chromosome

    def __len__(self) -> int:
        return self.__ALLELES
    
    def get_allele_size(self) -> int:
        return self.__DIM

    def mutate(self, count:int=0, percepts=[1, 1, 1, 1]):
        for i in range(self.__ALLELES):
            for _ in range(count):
                self.__chromosome[i].invert(randint(0, self.__DIM-1))

    def cross(self, other: 'Individual', mutate:float=None) -> 'Individual':
        if len(self) != len(other): raise NOT_EQUAL

        if self.__DIM != len(other.get_allele()[0]) : raise NOT_EQUAL

        # partner = other.get_allele()
        child = list(map( lambda x: Allele(self.__DIM, data=''.join(x)), self.__GETCHILDCHROMO( self, other )))

        child = np.array(child)
        child = Individual(self.__ALLELES, self.__DIM, child, __GETCHILDCHROMO=self.__GETCHILDCHROMO, MODULATOR=self.__MODULATOR)
        if mutate is not None: child.mutate(count=int(self.__DIM * mutate) )

        return child

    def __init__(self, __ALLELES: int, __DIM: int, __chromosome=None, mutate=None, __GETCHILDCHROMO:Callable=RANDOMCHROMO, MODULATOR: Callable=None, **kwargs):
        if __ALLELES < 2: raise SIZE_ERROR(__ALLELES, "ALLELES")
        if __DIM < 2  : raise SIZE_ERROR(__DIM, "DIM")
        if __chromosome is not None and (len(__chromosome) != __ALLELES or \
            len(__chromosome[0]) != __DIM): raise NOT_EQUAL(__ALLELES)

        self.__GETCHILDCHROMO = __GETCHILDCHROMO
        self.__MODULATOR = MODULATOR

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
    __MAXES: Tuple[float] # max_min if limits are equal
    __population: List[Individual] # List of Individuals
    __function: Callable
    __gen_res: List[Tuple[float]]
    __best: Tuple[Tuple, float]
    
    __LIMITS: Tuple[Tuple[float]] # table of maximum values for every dimension

    def get_accuracy(self):
        return self.__ACCURACY
    
    def get_best(self):
        return self.__best

    def __init__(self, populus_size: int, dim: int, dels: int, function: Callable, maxes: Tuple[float]=None, limits:Tuple[Tuple[float]]=None, IndividualArgs:dict={}, MODULATOR:Callable=None, KEYFUNCTION: Callable=None,reverse=True, **kwargs):
        
        self.__ACCURACY = math.ceil(math.log2(dels))# + 1
        
        self.__best = None
        # self.__mod = (max(maxes) - min(maxes)) / dels
        self.__DIM = dim
        
        self.IndividualArgs = IndividualArgs
        self.MODULATOR = MODULATOR

        self.__LIMITS = limits
        self.__MAXES = maxes

        self.reverse = reverse
        self.KEYFUNCTION = KEYFUNCTION

        self.__function = function
        self.__gen_res = None
        # print('MAXES', self.__MAXES)

        if self.__LIMITS is None:
            self.step = 1.0 * ( max(self.__MAXES) - min(self.__MAXES) ) / ((2**(self.__ACCURACY) - 1) + 1)
        else:
            self.step = []
            for i in range(len(self.__LIMITS)):
                self.step.append(
                    1.0 * ( max(self.__LIMITS[i]) - min(self.__LIMITS[i]) ) / ((2**(self.__ACCURACY) - 1) + 1)
                )

        self.__population = []
        for _ in range(populus_size):
            self.__population.append(
                Individual(
                    self.__DIM,
                    self.__ACCURACY,
                    MODULATOR=self.MODULATOR,
                    mutate=0.7,
                    **self.IndividualArgs
                )
            )
    
    def __modulate(self, x: int, index:int=0) -> 'np.ndarray[Allele]':
        # print("STEP >", f"[ {x} ] ", 2**self.__ACCURACY - 1, x*step)
        if self.__LIMITS is None:
            return self.MIN + (x * self.step)
        else:
            return self.__LIMITS[index][0] + (x * self.step[index])
        
    def count(self) -> list:
        #if self.__gen_res is None:
            #print("AAAAAAA")
        self.__best = None #(tuple(), self.__gen_res[0][0])

        self.__gen_res = []
        for i in range(len(self.__population)):
            #print(self.__population[i].get_vector())
            populus = self.__population[i].get_vector()
            if self.MODULATOR is None:
                xs = np.array([ self.__modulate(populus[x], index=x) for x in range(len(populus)) ])
            else:
                xs = np.array(self.MODULATOR(xs))
            #print(xs)
            z = self.__function(xs)
            self.__gen_res. append(
                (*xs, z) # (1, 2, 3, 4, 12412)
            )
            # if self.__best is None or (self.__best[1] >= z and not self.reverse) or (self.__best[1] <= z and self.reverse):
            #     self.__best = (tuple(xs), z)

        self.__best = sorted(
            enumerate(self.__population),
            key=lambda x: self.__gen_res[x[0]][-1] + self.KEYFUNCTION(self.__gen_res[x[0]][:-1]),
            reverse=self.reverse
        )[0]

        self.__best = (tuple(
            [
                self.__modulate(i)
                for i in self.__best[1].get_list()
            ]
        ), self.__gen_res[self.__best[0]])

        # print(self.__best)
        return [] + self.__gen_res

    def next(self):
        s = sorted(
            enumerate(self.__population),
            key=lambda x: self.__gen_res[x[0]][-1] + self.KEYFUNCTION(self.__gen_res[x[0]][:-1]),
            reverse=self.reverse
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