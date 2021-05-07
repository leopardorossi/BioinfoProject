from .input_parsing import read_kmers_countings
from .counting import Counter

if __name__ == '__main__':
    
    kmerData = read_kmers_countings("python_approach/resources/out.txt")
    kmerData = Counter.orderKmers(2, kmerData)

    