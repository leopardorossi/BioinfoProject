import operator

from enum import Enum
from collections import OrderedDict

class KmerData:
    def __init__(self, kmer, freq, family):
        self.kmer = kmer
        self.freq = freq
        self.family = family
    
    def __str__(self):
        return f"Kmer: {self.kmer}\tFreq: {self.freq}\tFamily: {self.family}"

class Counter:
    
    class SortingMethods(Enum):
        LEXICOGRAPHICAL=1
        FAMILY=2

    def __init__(self):
        pass

    """
        Orders the given array that contains data about kmers countings
        @param data The array that has to be stored
    """
    @staticmethod
    def orderKmers(method, data):
        if method == Counter.SortingMethods.LEXICOGRAPHICAL:
            return sorted(data, key=lambda kData: kData.kmer)
        elif method == Counter.SortingMethods.FAMILY:
            return sorted(data, key=lambda kData: kData.family)

    @staticmethod    
    def mapFirstLastBases(firstBase, lastBase):
        basesMap = {
            "AA": 1, 
            "AC": 2, 
            "AG": 3, 
            "AT": 4, 
            "CA": 5, 
            "CC": 6, 
            "CG": 7, 
            "CT": 8,
            "TA": 9,
            "TC": 10,
            "TT": 11,
            "TG": 12,
            "GA": 13,
            "GC": 14,
            "GT": 15,
            "GG": 16
        }
        key = firstBase + lastBase
        return basesMap[key]