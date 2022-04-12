
from .errors import *

class Allele:

    __data: str # Str of 0 and 1 - array of bits of gene.\
                # First bit describes a sign of <integer> interpretation.
                # Max value is (2^__DIM - 1)
    __DIM: int  # Fixed size of Allele (Genes count)

    def set(self, data):
        self.__data = data

    def __len__(self) -> int:
        return self.__DIM

    def __str__(self) -> str:
        return f"<Allele: {self.__DIM}> ({self.__data});"

    def __repr__(self) -> str:
        return self.__str__()
    
    def get_dumpstr(self) ->str:
        return "" + self.__data

    def get_data(self) -> int:
        return int(self.__data, 2)# * (-1 if self.__data[0] == "1" else 1)

    def get_bit(self, n: int) -> int:
        if self.__DIM <= n: raise INDEX_ERR
        return self.__data[n]

    def swap(self, r: int, l: int):
        if self.__DIM <= r or self.__DIM <= l or l < 0 or r < 0: raise INDEX_ERR

        d = str(list(self.__data))
        d[l], d[r] = d[r], d[l]
        self.__data = ''.join(d)
    
    def invert(self, n: int):
        # Inverts n bit from gen
        if self.__DIM <= n or n < 0: raise INDEX_ERR
        self.__data = self.__data[:n] + \
            str((int(self.__data[n]) + 1) % 2) + \
            (self.__data[n+1:] if (n+1) < self.__DIM else "")

    def insert(self, n: int, value: int):
        # Insets value to n bit of gen
        if value > 1 or value < 0: raise ValueError("Bit can be 0 or 1 only")
        self.__data = self.__data[:n] + str(int(value)) + \
            (self.__data[n+1:] if (n+1) < self.__DIM else "")

    def __init__(self, size:int, data=None):
        if size < 2:
            raise SIZE_ERROR(size, "size")

        self.__DIM = size

        if data is None:
            self.__data = '1' * (self.__DIM)
        else:
            self.__data = data
    
    def clear(self):
        self.__data = '1' * (self.__DIM)