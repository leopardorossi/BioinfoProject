from .counting import Preprocessing, KmerData

from functools import reduce


def read_kmers_counting(input_file):
    # Define the array that will contain the counting
    counting = []
    # Define an array that contains the sizes of each family
    families_sizes = []
    for i in range(0, 16):
        families_sizes.append(0)

    # Open the kmer counting file
    with open(input_file, "r") as f:
        # Loop over the lines to retrieve the counting
        for line in f:
            # Split the current line in order to retrieve the kmer and its counting
            kmer, freq = line.split("\t")
            # Produce an integer number that identifies the family of the kmer
            family = Preprocessing.map_first_last_bases(kmer[0], kmer[len(kmer) - 1])
            # Increment the size of the family
            families_sizes[family-1] += 1
            # Build an item to store all the information together
            freq = int(freq.rstrip("\n"))
            counting.append(KmerData(kmer, freq, family))

    return counting, families_sizes


def generate_seed_from_pattern(pattern, kmer_length):
    # Split the pattern into the different groups
    groups = pattern.split("-")
    pattern_length = reduce(lambda a, b: int(a)+int(b), groups)

    if pattern_length > kmer_length:
        raise Exception("The pattern length must be less or equal to the kmer length")

    # Generate the seed
    pattern = ""
    for i, group in enumerate(groups):
        if (i+1) % 2 == 0:
            zeros = "0" * int(group)
            pattern += zeros
        else:
            ones = "1" * int(group)
            pattern += ones

    # Check if the seed is palindrome
    reverse = pattern[::-1]
    if pattern != reverse:
        raise Exception("The seed must be palindrome")

    return pattern




