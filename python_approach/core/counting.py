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
        LEXICOGRAPHICAL = 1
        FAMILY = 2

    @staticmethod
    def order_kmers(method, data):
        if method == Preprocessing.SortingMethods.LEXICOGRAPHICAL:
            return sorted(data, key=lambda k_data: k_data.kmer)
        elif method == Preprocessing.SortingMethods.FAMILY:
            return sorted(data, key=lambda k_data: k_data.family)

    @staticmethod    
    def map_first_last_bases(first_base, last_base):
        bases_map = {
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
        key = first_base + last_base
        return bases_map[key]


class CountingStrategy(metaclass=abc.ABCMeta):

    def __init__(self, seed, data):
        self.seed = seed
        self.data = data

    def execute(self):
        space_seed_result = self.apply_space_seed()
        self.create_hash_table(space_seed_result)
        self.count()

    def apply_space_seed(self):
        # Get the indexes of the ones inside the seed
        ones_positions = [i for i, c in enumerate(self.seed) if c == "1"]
        
        # Define the list that will contain the result of seed applying
        result = []
        # Loop over the data and apply the seed to the kmer
        for kData in self.data:
            mask_applied_result = ""
            for index in ones_positions:
                mask_applied_result += kData.kmer[index]
            
            result.append(KmerData(mask_applied_result, kData.freq, kData.family))
        
        return result

    @abc.abstractmethod
    def create_hash_table(self, data):
        pass

    @abc.abstractmethod
    def count(self):
        pass


class LexicoGraphicalCountingStrategy(CountingStrategy):

    def __init__(self, seed, data):
        super().__init__(seed, data)

    def create_hash_table(self, data):
        pass

    def count(self):
        pass


class FamilyCountingStrategy(CountingStrategy):

    def __init__(self):
        pass

    def apply_space_seed(self):
        pass

    def create_hash_table(self, data):
        pass

    def count(self):
        pass
