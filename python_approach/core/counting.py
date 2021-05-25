import csv
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

    def __init__(self, seed, data, families_sizes, output_path):
        self.seed = seed
        self.data = data
        self.families_sizes = families_sizes
        self.output_path = output_path

    def execute(self):
        space_seed_result = self.apply_space_seed()
        hash_table = self.create_hash_table(space_seed_result)
        counting = self.count(hash_table)
        self.store_output_file(counting)
        self.generate_stats_file()

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

    def generate_stats_file(self):
        stats = {}
        with open(f"{self.output_path}/result_{self.description()}.txt", "r") as f:
            for line in f:
                kmer, freq = line.split("\t")
                freq = int(freq)
                if freq in stats.keys():
                    stats[freq] += 1
                else:
                    stats[freq] = 1

        with open(f"{self.output_path}/profile_{self.description()}.csv", "w") as f:
            f_writer = csv.writer(f, delimiter=";")
            f_writer.writerow(["Frequency", "Count"])
            for key in stats.keys():
                f_writer.writerow([key, stats[key]])

    def store_output_file(self, counting):
        with open(f"{self.output_path}/result_{self.description()}.txt", "w") as f:
            for key in counting.keys():
                f.write(f"{key}\t{counting[key]}\n")


class LexicoGraphicalCountingStrategy(CountingStrategy):

    def __init__(self, seed, data, output_path):
        super().__init__(seed, data, [], output_path)

    def create_hash_table(self, data):
        # With this method we do not have to create a real hash table, but
        # we only have to order the masked kmers to count them
        return Preprocessing.order_kmers(Preprocessing.SortingMethods.LEXICOGRAPHICAL, data)

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

    def create_hash_table(self, data):
        hash_table = {}

        start = 0
        for i, size in enumerate(self.families_sizes):
            if size > 0:
                hash_table[i+1] = data[start:(start+size)]
                start += size

        return hash_table

    def description(self):
        return f"{self.seed}_family"

