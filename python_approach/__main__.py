from python_approach.core import *

if __name__ == '__main__':

    kmer_data, families_sizes = read_kmers_counting("resources/out.txt")
    kmer_data = Preprocessing.order_kmers(Preprocessing.SortingMethods.FAMILY, kmer_data)

    # First seed: 1001001001001001001001001001
    # Second seed: 1000000000000000000000000001
    cStrategy = FamilyCountingStrategy("1001001001001001001001001001", kmer_data, families_sizes)
    cStrategy.execute()

    generate_stats_file("result.txt", "stats.txt")
