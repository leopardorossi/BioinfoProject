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

    def __init__(self, seed, data, families_sizes):
        self.seed = seed
        self.data = data
        self.families_sizes = families_sizes

    def execute(self):
        space_seed_result = self.apply_space_seed()
        hash_table = self.create_hash_table(space_seed_result)
        counting = self.count(hash_table)
        self.store_output_file(counting)

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

    def count(self, table):
        counting = {}

        for key in table.keys():
            for k_data in table[key]:
                # Check if the current kmer already exists in the hash table
                if k_data.kmer in counting.keys():
                    counting[k_data.kmer] += k_data.freq
                else:
                    counting[k_data.kmer] = k_data.freq

        return counting

    def store_output_file(self, counting):
        with open("resources/result.txt", "w") as f:
            for key in counting.keys():
                f.write(f"{key}\t{counting[key]}\n")


class LexicoGraphicalCountingStrategy(CountingStrategy):

    def __init__(self, seed, data):
        super().__init__(seed, data, [])

    def create_hash_table(self, data):
        pass


class FamilyCountingStrategy(CountingStrategy):

    def __init__(self, seed, data, families_sizes):
        super().__init__(seed, data, families_sizes)

    def create_hash_table(self, data):
        hash_table = {}

        start = 0
        for i, size in enumerate(self.families_sizes):
            if size > 0:
                hash_table[i+1] = data[start:(start+size)]
                start += size

        return hash_table


def generate_stats_file(count_file, out_file):
    stats = {}
    with open(f"resources/{count_file}", "r") as f:
        for line in f:
            kmer, freq = line.split("\t")
            if freq in stats.keys():
                stats[freq] += 1
            else:
                stats[freq] = 1

    with open(f"resources/{out_file}", "w") as f:
        for key in stats.keys():
            f.write(f"Freq: {key} Count:{stats[key]}\n")

