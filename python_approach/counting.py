import operator
import abc

from enum import Enum
from collections import OrderedDict

class KmerData:
    def __init__(self, kmer, freq, family):
        self.kmer = kmer
        self.freq = freq
        self.family = family
    
    def __str__(self):
        return f"Kmer: {self.kmer}\tFreq: {self.freq}\tFamily: {self.family}"

class Preprocessing:

    class SortingMethods(Enum):
        LEXICOGRAPHICAL=1
        FAMILY=2

    @staticmethod
    def orderKmers(method, data):
        if method == Preprocessing.SortingMethods.LEXICOGRAPHICAL:
            return sorted(data, key=lambda kData: kData.kmer)
        elif method == Preprocessing.SortingMethods.FAMILY:
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

class CountingStrategy(metaclass=abc.ABCMeta):

    def __init__(self, seed, data):
        self.seed = seed
        self.data = data

    def execute(self):
        spaceSeedResult = self.applySpaceSeed()
        self.createHashTable(spaceSeedResult)
        self.count()

    def applySpaceSeed(self):
        print(self.seed)
        # Get the indexes of the ones inside the seed
        onesPositions = [i for i, c in enumerate(self.seed) if c == "1"]
        
        # Define the list that will contain the result of seed applying
        result = []
        # Loop over the data and apply the seed to the kmer
        for kData in self.data:
            maskAppliedResult = ""
            for index in onesPositions:
                maskAppliedResult += kData.kmer[index]
            
            result.append(KmerData(maskAppliedResult, kData.freq, kData.family))
        
        return result

    @abc.abstractmethod
    def createHashTable(self, data):
        pass

    @abc.abstractmethod
    def count(self):
        pass

class LexigoGraphicalCountingStrategy(CountingStrategy):

    def __init__(self, seed, data):
        super().__init__(seed, data)

    def createHashTable(self, data):
        pass

    def count(self):
        pass

class FamilyCountingStrategy(CountingStrategy):

    def __init__(self):
        pass

    def applySpaceSeed(self):
        pass

    def createHashTable(self, data):
        pass

    def count(self):
        pass