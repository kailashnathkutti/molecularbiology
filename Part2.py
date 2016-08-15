'''
Created on 30.07.2016

@author: kuttik
'''
import random
from _collections import OrderedDict
# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    innerMotif = Count(Motifs)
    keys = innerMotif.keys()
    profile = {}
    count =0;

    for key in keys:
        profile[key]=[]
    for (key,value) in innerMotif.items():
        for i in value:
            count=count+1
            profile[key].append(i/t)
    # insert your code here
    #print profile
    return profile
#print (Profile(["ACGTTA",
#                 "AGGTGA",
#                 "ACGTTT",
#                 "ATGCTA",
#                 "TCGCGA",
#                 "AAGCTA"]))
#print (Profile(["GTC","CCC","ATA","GCT"]))

def Consensus(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
#print(Consensus(["ACGTTA"]))

#Profile(["AACGTA","CCCGTT"])


def Score(Motifs):
    z = zip(*Motifs)
    totalscore = 0
    for string in z:
        score = len(string)-max([string.count('A'),string.count('C'), string.count('G'), string.count('T')])
        totalscore += score
    return totalscore    
#print (Consensus(["ACGTTA","AAGAGA","AGGTGA","AGGTCA","ACGCGA","ATGCTA"]))
# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    chars = set('ACTG')
    # insert your code here
    retProb=1
    textPos=0
    for profileChar in Profile.keys():
        if profileChar not in chars:
            return 0.0            
    for symbol in Text:
        profileRow = Profile[symbol];
        retProb=retProb*profileRow[textPos]
        textPos=textPos+1
    return retProb

#Profile1 = { 'A': [0.4,0.3,0.0,0.1,0.0,0.9], 'C': [0.2,0.3,0.0,0.4,0.0,0.1],'G': [0.1,0.3,1.0,0.1,0.5,0.0],'T': [0.3,0.1,0.0,0.4,0.5,0.0] }
#Text = ["ACGTTA","AGGTGA","ACGTTT","ATGCTA","TCGCGA","AAGCTA"]
#Text = "CAGTGA"
#print(Pr(Text,Profile1))

def ProfileMostProbablePattern(Text, k, Profile):
    l = len(Text)
 
    pr_most = 0.0
    mostProbablePattern = Text[0:k]
 
    for j in range(l-k+1):
        k_mer = Text[j:j+k]
        #print (k_mer)
        pr = Pr(k_mer, Profile)
        if pr > pr_most:
            pr_most = pr
            mostProbablePattern = k_mer
 
    return mostProbablePattern 
# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)


def GreedyMotifSearch(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
 
    # length of the first Dna string
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))      
 
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
 
    return BestMotifs

# Copy the ten strings occurring in the hyperlinked DosR dataset below.
#Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
#Dna = ["GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG"]
#Dna = ["GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG"]
#Dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
# set t equal to the number of strings in Dna and k equal to 15
#Dna = ["GGCGTTCAGGCA"]
#k=3
#t=len(Dna)
# Call GreedyMotifSearch(Dna, k, t) and store the output in a variable called Motifs
#Motifs = GreedyMotifSearch (Dna,k,t)
# Print the Motifs variable
#print (Motifs)
#print (min(enumerate(Motifs[0])))
#print (min(range(len(Motifs[0]), key=Motifs[0].__getitem__)))


# Motifs.sort()
# #print(Motifs)
# print(Motifs[0][1])
# 
# print (Score(Motifs[0][1]))

#print(Motifs)
#VarMotifs.sort()
#print(Motifs)
#print(Motifs)
#print (Motifs)
# Print the Motifs variable
#print (Score(Motifs))

def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} # initializing the count dictionary
    # your code here
    
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    #print(count)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    #print(count)            
    return count
#print (CountWithPseudocounts(["AACGTA","CCCGTT","CACCTT","GGATTA","TTCCGG"]))

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    #Number of observations is increased by 1 magnitude as per la place rule
    t = len(Motifs)+4
    k = len(Motifs[0])
    profile = {} # output variable
    innerMotif = CountWithPseudocounts(Motifs)
    keys = innerMotif.keys()
    #print(innerMotif, keys)
    for key in keys:
        profile[key]=[]
    for (key,value) in innerMotif.items():        
        for i in value:
            #print(key,value,i,t,float(i/t))
            profile[key].append(float(i/t))
    # insert your code here
    #print profile
    # your code here
    return profile

#print (ProfileWithPseudocounts(["AACGTA","CCCGTT","CACCTT","GGATTA","TTCCGG"]))

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
 
    # length of the first Dna string
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))      
 
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
 
    return BestMotifs
#k = 15
#t = 10
#Dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]

#print(len(Dna))
#print(GreedyMotifSearchWithPseudocounts(Dna,k,t))
#print (Score(GreedyMotifSearchWithPseudocounts(Dna,k,t)))
#def Motifs(Profile, Dna):
    # insert your code here


#print(ProfileMostProbablePattern(Text,4,Profile1))


def Motifs(Dna,Profile):   
    t = len(Dna[0])
    k=len(Profile['A'])    
    n = len(Dna)    
    BestMotifs = {}
    retMotifs =[]
    for i in range(n):        
        BestMotifs[i]=[] 
        t = len(Dna[i])
        maxPr=0.0
        for j in range(0, t):            #P = Profile(Motifs[0:j])
            dnaSubstr=Dna[i][j:j+k]            
            if len(dnaSubstr) >=k :                
                if (Pr(dnaSubstr,Profile) > maxPr):
                    maxPr=Pr(dnaSubstr,Profile);
                    BestMotifs[i]=[dnaSubstr,maxPr]
                    
    for items in BestMotifs.items():        
        retMotifs.append(items[1][0])
        
           
    return retMotifs
    
#Test
Dna = ["TGACGTTC","TAAGAGTT","GGACGAAA","CTGTTCGC"]
print(Motifs(Dna,Profile(["TGA","GTT","GAA","TGT"])))
#Profile1 = { 'A':[0.8,0.0,0.0,0.2], 'C': [0.0,0.6,0.2,0.0],'G': [0.2,0.2,0.8,0.0],'T': [0.0,0.2,0.0,0.8] }
#Dna = ["TTACCTTAAC","GATGTCTGTC","ACGGCGTTAG","CCCTAACGAG","CGTCAGAGGT"]
#Dna = ["GATGTCTGTC","CGTCAGAGGT"]
#print(Motifs(Dna,Profile1))

def RandomMotifs(Dna, k, t):
    n = len(Dna)
    retMotifs =[]
    for i in range(n):
        randStart=int(random.randint(1,3))        
        if randStart > 0 and randStart <= k:
            retMotifs.append(Dna[i][randStart:randStart+k])
    return retMotifs

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)        
        M = Motifs(Dna,Profile)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs


### DO NOT MODIFY THE CODE BELOW THIS LINE ###
def RepeatedRandomizedMotifSearch(Dna, k, t,runs):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(runs):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs

k = 3
t = 4
runs=1
#Dna = ["ATGAGGTC","GCCCTAGA","AAATAGAT","TTGTGCTA"]
#print(RandomizedMotifSearch(Dna,3,4))
#BestMotifs = RepeatedRandomizedMotifSearch(Dna,k,t,runs)
#print (BestMotifs)
# print (Score(BestMotifs))




def Normalize(Probabilities):
    sum=0
    normList={}
    for value in Probabilities.items():
        sum =sum+value[1]
    for (key,value) in Probabilities.items():
        normList[key]=value/sum
        print(value/sum)
    return normList

#Probabilities = {'AA': 0.22, 'TT': 0.54, 'CC': 0.58, 'GG': 0.36,  'AT': 0.3 }
#print('Normalized', Normalize(Probabilities))

def WeightedDie(Probabilities):
    kmer = '' # output variable
    # your code here    
    kmerProbabilityRangeMap = {}
    kmerProbSum=sum(Probabilities.values())
    kmerProbInterimSum=0.0
    counter=0    
    probList = list(Probabilities.values());    
    probListLen=len(probList)-1
    for key,value in Probabilities.items():                
        if counter == 0:
            kmerProbabilityRangeMap[key]=[0,float(value)]
            counter=counter+1
            kmerProbInterimSum=float(value)
        elif counter > 0 and counter < probListLen:
            tmp=kmerProbInterimSum+probList[counter]
            kmerProbabilityRangeMap[key]=[kmerProbInterimSum,tmp]
            kmerProbInterimSum=kmerProbInterimSum+probList[counter]
            counter=counter+1
        elif counter == probListLen:
            kmerProbabilityRangeMap[key]=[kmerProbInterimSum,kmerProbSum]
     
    while True:
    #for lu in range(0,10000,1):
        randomRoundPr = random.uniform(0, 1)
        for key,value in kmerProbabilityRangeMap.items():            
            if randomRoundPr >= value[0] and randomRoundPr <=value[1]:
                #print(key,randomRoundPr)
                kmer =key
                #break
        break                
    return kmer


#print (WeightedDie(Probabilities))


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

#Text="AAACCCAAACCC"
#profile={'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
#k=2
#print(ProfileGeneratedString(Text, profile, k))

def getIndexOfChar(chr):
    return "ACGT".index(chr)

def getRandomMotifs( dna, k ):
    return [createRandomKmer(seq,k) for seq in dna]

def createRandomKmer( seq, k ):
    start = random.randint(0, len(seq)-k)
    return seq[start:start+k] 

def Score( motifs ):
    score = 0
    for count in countsFromMotifs(motifs,0):
        score += sum(count) - max(count)
    return score

def countsFromMotifs( motifs, initCount ):
    k = len(motifs[0])
    for motif in motifs:
        assert k == len(motif)
    counts = []
    for i in range(k):
        currCount = [initCount] * 4
        for motif in motifs:
            currCount[getIndexOfChar(motif[i])] += 1
        counts.append(currCount)
    return counts


def Random( pr ):
    total = float(sum(pr))
    r = random.random()
    partialSum = 0.0
    for i in range(len(pr)):
        partialSum += pr[i]
        if partialSum/total >= r:
            return i
    return -1

def findProfile( motifs, initCount ):    
    counts = countsFromMotifs(motifs,initCount)
    profile = []
    for count in counts:
        total = float(sum(count))
        probs = [c/total for c in count]
        profile.append(probs)
    return profile

def findProbKmerInProfile( kmer, profile ):
    prob = 1.0
    for i in range(len(kmer)):
        prob *= profile[i][getIndexOfChar(kmer[i])]
    return prob

def GibbsSampler( Dna, k,t, N ):
    BestMotifs = []
    motifs = getRandomMotifs(Dna,k)
    BestMotifs = motifs
    bestScore = Score(BestMotifs)
    for step in range(N):
        i = random.randint(0, t-1)
        profile = findProfile(motifs[:i] + motifs[i+1:], 1)
        pr = [findProbKmerInProfile(Dna[i][s:s+k],profile) for s in range(len(Dna[i])-k+1)]
        s = Random(pr)
        motifs[i] = Dna[i][s:s+k]
        score = Score(motifs)
        if score < bestScore:
            BestMotifs = motifs[:]
            bestScore = score
    
    return BestMotifs



#Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
# Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC"
# ,"CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG"
# ,"ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC"
# ,"GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC"
# ,"GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG"
# ,"CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA"
# ,"GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA"
# ,"GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG"
# ,"GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG"
# ,"TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
# 
# k=15
# t=10
# N=100
# print(GibbsSampler(Dna, k, t, N))
# print(Score(GibbsSampler(Dna, k, t, N)))