# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates

# Input:  A list Items
# Output: A list containing all objects from Items without duplicates
def remove_duplicates(listing):
    ItemsNoDuplicates = [] # output variable
    n = []    
    # print listing
    for part in range(len(listing)):
        exists = 0
        for part2 in range(len(n)):
            if listing[part] == n[part2]:
                exists = 1
        if exists == 0:
            n.append(listing[part])  
            ItemsNoDuplicates = n
    return ItemsNoDuplicates

# Input:  A string Text and an integer k
# Output: CountDict(Text, k)
# HINT:   This code should be identical to when you last implemented CountDict
def CountDict(Text, k):
    Count = {} # output variable
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count

# Input:  Strings Pattern and Text
# Output: The number of times Pattern appears in Text
# HINT:   This code should be identical to when you last implemented PatternCount
def PatternCount(Pattern, Text):
    count = 0 # output variable
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def complement(Nucleotide):
    length = len(Nucleotide)    
    if length < 0 :
        output[0] = '0'
    elif length > 0:
        output=[]
        outputStr=''
        comp = ''
        counter = length -1        
        for c in range(len(Nucleotide)) :            
            output.append(Nucleotide[counter])
            outputStr=outputStr+(Nucleotide[counter])
            counter = counter-1            
    comp = outputStr
    return outputStr


def ReverseComplement(Pattern):
    original=complement(Pattern)
    length = len(Pattern)    
    if length < 0 :
        outputStr=''
    elif length > 0:
        outputStr=''
        revComp = ''
        for c in range(len(Pattern)) :                        
            if original[c] == 'A' :
                outputStr=outputStr+'T'
            elif original[c] == 'T' :
                outputStr=outputStr+'A'
            elif original[c] == 'G' :
                outputStr=outputStr+'C'
            elif original[c] == 'C' :
                outputStr=outputStr+'G'
    revComp = outputStr             
    return revComp

### DO NOT MODIFY THE CODE BELOW THIS LINE ###
import sys

#print(' '.join(FrequentWords('ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC',10)))


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def Skew(Genome):
    skew = {} #initializing the dictionary
    skewList = [0]    
    for i in range(0, len(Genome)):
        if Genome[i] == 'C':
            skewList.append(skewList[i] - 1)
        elif Genome[i] == 'G':
            skewList.append(skewList[i] + 1)
        else:
            skewList.append(skewList[i])
    for idx, val in enumerate(skewList):
        skew[idx]=val
        #print(idx, val)
    
    return skew

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions =[]
    skew = Genome
    s = Skew(Genome)
    minimumValue = min(s.values())

    for (k,v) in s.items():
        if v == minimumValue:
            positions.append(k)
       # your code here
    return positions

skew = Skew ('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
#print (MinimumSkew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'))
#print (MinimumSkew('GATACACTTCCCGAGTAGGTACTG'))
#print(' '.join([str(skew[i]) for i in sorted(skew.keys())]))
#print (Skew('GATACACTTCCCGAGTAGGTACTG'))


def HammingDistance(p, q):
    count=0
    if len(p) == len(q):
        for c in range(0,len(p)):
            if p[c] != q[c]:
                count=count+1
    return count
#print (HammingDistance('CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA','CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG'))            

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)],Pattern) <= d:
            positions.append(i)
    
    return positions
#print ApproximatePatternMatching('ATTCTGGA','CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT',3)
print (ApproximatePatternMatching('ATCG','ATCG',4))


