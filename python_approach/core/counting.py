import abc
import tracemalloc

from enum import Enum


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
            data.sort(key=lambda k_data: k_data.kmer)
        elif method == Preprocessing.SortingMethods.FAMILY:
            data.sort(key=lambda k_data: k_data.family)

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

    def __init__(self, seed, data, families_sizes, output_path):
        self.seed = seed
        self.data = data
        self.families_sizes = families_sizes
        self.output_path = output_path

    def execute(self):
        # tracemalloc.start()

        print("Hash table definition...")
        hash_table = self.create_hash_table()
        # current, peak = tracemalloc.get_traced_memory()
        # print(f"Current memory usage is {current / 10 ** 6}MB; Peak was {peak / 10 ** 6}MB\n")

        print("Counting...")
        counting = self.count(hash_table)
        # current, peak = tracemalloc.get_traced_memory()
        # print(f"Current memory usage is {current / 10 ** 6}MB; Peak was {peak / 10 ** 6}MB\n")

        # tracemalloc.stop()

        return counting

    @abc.abstractmethod
    def create_hash_table(self):
        pass

    @abc.abstractmethod
    def description(self):
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


class LexicoGraphicalCountingStrategy(CountingStrategy):

    def __init__(self, seed, data, output_path):
        super().__init__(seed, data, [], output_path)

    def create_hash_table(self):
        # With this method we do not have to create a real hash table, but
        # we only have to order the masked kmers to count them
        return self.data

    def count(self, table):
        # Thanks to the ordering we know that all equal masked kmer are close together, so we
        # exploit this while counting
        counting = {}

        for item in table:
            if item.kmer in counting.keys():
                counting[item.kmer] += item.freq
            else:
                counting[item.kmer] = item.freq

        return counting

    def description(self):
        return f"{self.seed}_lexi"


class FamilyCountingStrategy(CountingStrategy):

    def __init__(self, seed, data, families_sizes, output_path):
        super().__init__(seed, data, families_sizes, output_path)

    def create_hash_table(self):
        hash_table = {}

        start = 0
        for i, size in enumerate(self.families_sizes):
            if size > 0:
                hash_table[i+1] = self.data[start:(start+size)]
                start += size

        return hash_table

    def description(self):
        return f"{self.seed}_family"

