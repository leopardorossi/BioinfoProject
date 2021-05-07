from .input_parsing import read_kmers_countings
from .counting import Preprocessing, LexigoGraphicalCountingStrategy

if __name__ == '__main__':
    
    kmerData = read_kmers_countings("python_approach/resources/out.txt")
    kmerData = Preprocessing.orderKmers(Preprocessing.SortingMethods.FAMILY, kmerData)

    cStrategy = LexigoGraphicalCountingStrategy("1001001001001001001001001001", kmerData)
    cStrategy.execute()
    