# %%
from typing import Iterator
from itertools import chain
from functools import partial, reduce
import itertools
from itertools import combinations
import numpy as np
import re
from typing import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

#%%
# First practice - count seq
dna_seq = Seq("GCCGGCTGGTAGCCACTGCTTTTTCCCGAGCACCACAAATTATCATTCTGTGATCATTGTAGAATTCTCAACGGAGGTTTTTCCGAACGAAACCCCGTCGGACGCGGGGGTACCGGTAAACTGCTCCGGCAATGGCTTAATAACCGTAAACTCAGCCGGTTTACTTAGAAGCCCTAGCTACCTAACACCCACCGCGCCAATAACATTCTCTAAAGATAGTCACCCCAGGTCCCTTACCTGCATAGCCCTACGTTAAACCAGTTTTCTCGCACATCGAACCAGGAGTTATGCCCGCGTAGACCTGTATAGTGTAGTTGTAGTCGCAACGGGGCTTAAGAGCAAATCAGATTCCTAATGCTGGTCACAACCCGATGCACCTGGGATAATCCTGCTTACACGTTGTTCCTCCGACGTTAGCGAGGCCACGAACATCAAGTTTATCGGAGATGCGAACCATACCTAGATTTCCCCAACTCGATTACACTCGGATGAAAATTTACCAACCCGAACGTGTAGTATAAGATCACTCGTTCTTCTATGTGTCTATATGTGTAAGTCGGATTAACATACAGCAAATCACCAAAGTCCGGATACCGGCAGCAAAAGCCACGCTGCATAACCAGATCGATACCCGCTACTAGTGTCAACGACCTCCGTACCCTGTACTTTGAGCTTGGTCGTCGTCGCCGTCGATTAAAGGTGCCTGGCGGATCAGAGCTGTTTGGAATTAGTACTCCGTCTGATACCAAAGTTTACCCTCCTCCAAAGGGAAGTGATTTCCCTAGGAAATTGTCGAACCCTCCAAGTATGTTGTTGGATTATAATACTAGGAACGTGTGTAAACTTGGCAGCCAACTGACCCTCAGTCGGCGCGATGATAGTACGTGTACATACCCCTGCGGG")
num_A = dna_seq.count("A")
num_C = dna_seq.count("C")
num_G = dna_seq.count("G")
num_T = dna_seq.count("T")

print(num_A, num_C, num_G, num_T)

################################

# Best solution
def cout(s):
    return s.count("A"), s.count("G"), s.count("C"), s.count("T")

cout(dna_seq)

# %%
#Second practice - convert DNA to RNA
DNA = "TTGATATAGGTTGGCTCGCGTCAAAAAACGTGTCGTCTCACCGGGGCATGCGTTTTGGCCTATTCTCCGGGAGCATTCGATTATCCTTCCCTGCTTGTGCCACAGAGGTCCCCACGGATACTCGGACTGGGGATAAGCGTAGTGCTATCTGTACAAGTTTTGAATACCCCTGCTGGCATATAACTGCAGGATGTGCAATTAGAATCTCTGCAGTCATGGGAATTCGGTGATGACACGGAATGTTTCTAAAAACCGTGTGCCCGATCTCCGAAACGGTAAAAAGGCGACGATTACCGTTTCATCATCGACGGAACAGATGCGGACTCCTGTGGCTGCAGAAGCACGCGTGGCCCGTCTTGCGGTATCCCTGCTGGAATGCACCAGATCCACTCGCTTAAACACCACGTAGACGCGCCTATTACATGGCAGGCCAGGTTCGCAATAGTACAGACCACTTAGCAACTCGAAAAGAATAACATTGTGGTCCACTCCTGCGCTTAATGGGGTTGGTGTGCTTGCTCAAACTCGTTCGGCCGCACGTTTGCGTTCCCGTAAATTGCGTAGCGCAGAACCTAAATAATCGGTGCTTGCATCGTGGTCGTGGCGTAATTTACATGAACGCTTAGAGAACAGATGCCCGGTTCGGCTATTCTGCGTTGTTCCTACAGACTTTTGCCGCTGGTGACCTAGAGTACACGCGAGTCGAGATCAATGGGAGGGCCCGGTGACGGCATTTTGTAGACTGTAAGGAATTAGCTATAAAGCGTTGGGGAACAGACATGTCGCGATTGCGCTCAATAACGGTGATACAGTGGGCAATATCCATTTGGTAATAGCCGGGCTTCGTTAATGGTAGCCTTCAAGGGCGTAGCTTCGTCAGGTAACTGATGTGGTTCGGATGTAAAAATAGTGCAAAATGATCCGGCACATCTGGCGATTTAATAGGAAGGTGCAGAGAT"
RNA = DNA.replace("T", "U")
print(RNA)

################################

# Best solution
s = input()
print(s.replace('T', 'U'))
# %%
# Third practice
dna = input()
def trans(s):
    dictionaryS = {"A":"T", "C":"G","G":"C", "T":"A"}
    transTable = s.maketrans(dictionaryS)
    final = s.translate(transTable)
    print(final[::-1])
trans(dna)
################################

# Best solution
print(dna[::-1].translate(str.maketrans('ACGT', 'TGCA')))
# %%
# Fouth practice - rabbit count (fabonacci)
def fb(n, k):
    if n == 1 or n == 2:
        return 1
    else:
        return fb(n-2, k) * k + fb(n-1, k)
num, kk = int(input()), int(input())

fb(num,kk)
################################

# Best solution
def fib(n, k):
    a, b = 1, 1
    for i in range(2, n):
        a, b = b, k*a + b
    return b
fib(num, kk)
# %%
# Fifth practice - GC content in FASTA
## With biopython
records = list(SeqIO.parse("sample.fasta", "fasta"))
num = len(records)

gc_lists = [GC(records[i].seq) for i in range(0, num)]
maxx = max(gc_lists)
index = gc_lists.index(maxx)
label = records[index].id
print(label, "\n", maxx)

################################
## Without biopython
def parseChunk_GC(s):
    # parse label and seq
    seq = s.split('\n', maxsplit=1)
    # calcualte GC content
    GC = (seq[1].count("G") + seq[1].count("C"))/(len(seq[1]) - seq[1].count('\n'))
    return seq[0], GC

with open("sample.fasta") as f:
    lines = f.read()
    allDNA = lines.split('>')

    index, max = 0, 0
    for i in range(1, len(allDNA)):
        curr = parseChunk_GC(allDNA[i])[1]
        if curr > max:
            max = curr
            index = i
    label = parseChunk_GC(allDNA[index])[0]
    print(label, max)

# Best solution
with open('sample.fasta', mode='r', encoding='utf+8') as file:
    re = ''.join(file.read().split('\n')).split('>')
    re.remove('')
    d = {}
    for i in re:
        x = i[13:].count('C') + i[13:].count('G')
        d[i[:13]] = x * 100 / (len(i[13:]))
    print(max(d, key=d.get), max(d.values()))
# %%
# Sixth practice - count point mutation
with open("code.txt") as file:
    mut = file.read().split('\n')
    count = sum([i != j for i, j in zip(mut[0], mut[1])])
    print(count)

################################
# Best solution
print(sum(map(str.__ne__, *open("code.txt").read().split())))
# %%
# Seventh practice
k, m, n = int(input()), int(input()), int(input())
t = k+m+n
equation = k/t + m/t * (4*k+3*m+2*n-3)/(4*t-4) + n/t * (2*k+m)/(2*t-2)
print(equation)

################################

# Best solution
def firstLaw(k, m, n):
   N = float(k+m+n)
   return(1 - 1/N/(N-1)*(n*(n-1) + n*m + m*(m-1)/4.))

firstLaw(k, m, n)
# %%
# Eighth practice - DNA translation
## With biopython
seq = input()
print(Seq.translate(seq)[:-1])

################################
## Without biopython
table = "UUU F      CUU L      AUU I      GUU V\
 UUC F      CUC L      AUC I      GUC V\
 UUA L      CUA L      AUA I      GUA V\
 UUG L      CUG L      AUG M      GUG V\
 UCU S      CCU P      ACU T      GCU A\
 UCC S      CCC P      ACC T      GCC A\
 UCA S      CCA P      ACA T      GCA A\
 UCG S      CCG P      ACG T      GCG A\
 UAU Y      CAU H      AAU N      GAU D\
 UAC Y      CAC H      AAC N      GAC D\
 UAA Stop   CAA Q      AAA K      GAA E\
 UAG Stop   CAG Q      AAG K      GAG E\
 UGU C      CGU R      AGU S      GGU G\
 UGC C      CGC R      AGC S      GGC G\
 UGA Stop   CGA R      AGA R      GGA G\
 UGG W      CGG R      AGG R      GGG G"

seq = input()
codes = table.rstrip().split()
genes, letter = [], []
[genes.append(codes[i]) if i % 2 == 0 else letter.append(codes[i])
 for i in range(0, len(codes))]
check_table = dict(zip(genes, letter))
pase_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
result = ''.join([dict.get(check_table, item) for item in pase_seq][:-1])
print(result)

# Best solution
string = """UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G"""

coded = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
decoded = ''

traL = string.split()
traDict = dict(zip(traL[0::2], traL[1::2]))

for i in range(0, len(coded)-3, 3):
    decoded += traDict[coded[i:i+3]]

print(decoded)
# %%
# Ninth practice - find motif
s, motif = str(input()).split()
print([i+1 for i in range(len(s)) if s.startswith(motif, i)])

################################
s1, s2 = open('rosalind_subs.txt').read().split('\r\n')

for i in range(len(s1)):
    if s1[i:].startswith(s2):
        print(i+1)
# Best solution

# %%
# Tenth practice
with open('sample.fasta') as file:
    dna = ''.join(file.read().split('\n')).split('>')[1::]
    dna = [re.sub("[a-z]+_[0-9]+", '', i)[1::] for i in dna]

    res = [list(sub) for sub in dna]
    res_invert = list(zip(*res))

    count_m = [(item.count("A"), item.count("C"), item.count("G"), item.count("T")) for item in res_invert]
    dicts = {"0": "A", "1": "C", "2": "G", "3": "T"}
    A = [str(item.index(max(item))) for item in count_m]
    matchs = [dicts.get(item) for item in A]
    result_m = (np.asarray(count_m)).T
    s = ""
    for i in range(len(matchs)):
        s += matchs[i]
    print(s)
    for i in range(4):
        sto = dicts.get(str(i)) + " :"
        print(sto, end =" "), print(*(list(result_m[i])), " ")

################################

# Best solution
# %%
# Eleventh practice
def fib(passed_month: int, max_age: int, mature_month: int = 2):
    current_generation = [1] + [0] * (max_age-1)
    for i in range(passed_month-1):
        current_generation = [
            sum(current_generation[mature_month-1:])] + current_generation[:-1]
    return(current_generation)

print(sum(fib(85,16)))
################################
# Best solution
## Using the recursive function to calculate
def getChilds(n, m):  # find child population Cn
    if (n < 1):
        return 0
    else:
        if (n == 1):
            return 1
        else:
            sumChild = 0
            for i in xrange(2, m + 1):
                sumChild += getChilds(n - i, m)
            return sumChild


def getPop(n, m):  # find current population Fn
    return getChilds(n, m) + getChilds(n + 1, m)


# %%
# Twieth practice - Find the adjacency list
# check suffix of a string with the prefix of another string
def matches(s, t):
    if s[-3:] == t[:3]:
        return 1

records = list(SeqIO.parse("store.fasta", "fasta"))
for pairs in combinations(records, 2):
    #print(pairs[0].seq[-3:], pairs[1].seq[:3])
    if matches(pairs[0].seq, pairs[1].seq) == 1:
        print(pairs[0].id, pairs[1].id)

################################

# %%
# Thirtenth practice
import numpy
number_c = input().split()
# number_c = "18855 19614 16897 18945 16056 16489".split()
prob = [1, 1, 1, 3/4, 1/2, 0]
k = 2
sums = [float(number_c[n]) * prob[n] * k for n in range(6)]
print(sum(sums))

################################

# Best solution
results = print(sum([a*int(b) for a, b in zip([2, 2, 2, 1.5, 1, 0], input().split())]))

# %%
# Fourtinth practice - find the longest common string in a DNA string
records = list(SeqIO.parse("store.fasta", "fasta"))
maxlen = len(records)
sequences = list([records[i].seq for i in range(maxlen)])

def ngram(seq: str, n: int) -> Iterator[str]:
    return (seq[i: i+n] for i in range(0, len(seq)-n+1))


def allngram(seq: str) -> set:
    lengths = range(len(seq))
    ngrams = map(partial(ngram, seq), lengths)
    return set(chain.from_iterable(ngrams))


# sequences = ["brownasdfoersjumps",
#              "foxsxzxasis12sa[[#brown",
#              "thissasbrownxc-34a@s;"]

seqs_ngrams = map(allngram, sequences)
intersection = reduce(set.intersection, seqs_ngrams)
longest = max(intersection, key=len)  # -> brown
print(longest)







################################

# Best solution


# %%
# Fourtinth practice - find the longest common string in a DNA string

################################

# Best solution

# %%
# Fourtinth practice - find the longest common string in a DNA string

################################

# Best solution
