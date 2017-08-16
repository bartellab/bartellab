# this module contains information for 21U-RNA motif scoring

import motifModule

# THE MOTIF INTERFACE
# motif['len']: the length of the motif in nucleotides
# motif['apply'](seq): a lambda that applies the motif definition to the candidate seq 
#     and returns a numerical assessment.  what the numbers represent, and what scope of
#     values they can encompass, may change from case to case.
# motif['string']: is a string that represents the motif somehow
# motif['name']: is a string with a name for the motif


# pssm_motifs_A will be a list of the three string features for the
# TWO_PSSM i want to make
pssm_motifs_A = []


# the counts below are taken from Ruby et al 2006, and include
# loci from 21U-RNA-rich portions of the genome for which reads
# were observed.  Note that neither a length of 21nt nor a 5p U
# were requirements for a locus to be included in the conts below.
# sequences for the large motif were centered around the core motif
# as described in Ruby et al. 2007.

largeMotif = motifModule.PSSM("""
background .34  .34    .16    .16
3-6-06   A       T       C       G
0      2792    2464    710     628    
1      3122    2044    765     663    
2      3678    1587    869     460    
3      4129    1612    465     388    
4      4262    1636    223     473    
5      3675    2077    324     518    
6      2996    2663    633     302    
7      2798    2940    441     415    
8      1792    3930    439     433    
9      1507    3143    1354    590
10     1221    3604    1031    738
11     4079    1100    741     674
12      926     369     5184    115   
13      871     4951    435     337   
14      290     355     84      5865  
15      57      6246    259     32    
16      68      6488    26      12    
17      105     6340    94      55    
18      45      530     5948    71    
19      5543    425     118     508
20      2581    1292    1693    1028
21      2514    1652    1123    1305
22      2825    2426    457     886
23      2281    3297    562     454
24      1292    1704    1045    2553  
25      977     4218    713     686   
26      636     4372    1149    437   
27      1264    2423    1520    1387  
28      3185    2731    183     495   
29      308     4562    1215    509   
30      4545    1603    61      385   
31      3178    3239    89      88    
32      4509    914     422     749   
33      3240    2554    280     520   
""")



# the distance distribution was determined as described in
# Ruby et al. 2006.

distance = motifModule.DISTANCE_PSSM("""
dist   25   24   23   22   21   20   19   18   17   16
blah    6   28   82  388 1225 1979 1108   73   32    1
""")



# the small motif counts were obtained as described for the large motif
# above, but with sequences centered around the 5p position of the
# observed read (position 0 below).

smallMotif = motifModule.PSSM("""
background .34  .34    .16    .16
index       A     T     C      G
-4        1784  3216   969    907 
-3         497  4025  2097    259 
-2        4604   187    99   1987
-1        2395  2723   808    951 
0           58  6551   111    158 
1         2243  1152  1376   2107
""")







# the number of nucleotides downstream of the 21U-RNA 5p nucleotide to include
downLen = 1 ## because the small motif includes the first TWO nucs of the 21U

# the number of nucleotides downstream of the 21U-RNA 5p nucleotide to include;
# must include the length of the upstream portion of the small motif, the length
# of the big motif, and the max allowed distance between the two.
upLen = largeMotif['len'] + distance['max len'] + 4 ## four is the number of upstream nucs
                                                    ## in the small motif


# specific for this module; applies the motif with knowledge of the 5' end of
# the RNA that will be generated.
def applyScore(seq):
    if len(seq)!=1+upLen+downLen:
        raise ValueError("length of sequence was incorrect (" + str(len(seq)) + \
                         " instead of " + str(1+upLen+downLen) + "):\n" + seq + "\n")
    smallScore = smallMotif['apply'](seq[0-smallMotif['len']:])
    max_score = None
    for n in range(distance['min len'], distance['max len']+1):
        start = len(seq)-n-smallMotif['len']-largeMotif['len']
        end = start + largeMotif['len']
        new_score = largeMotif['apply'](seq[start:end]) + distance['apply'](n) + smallScore
        if max_score==None or new_score>max_score: max_score = new_score
    return max_score


# scores with the large motif, then finds the best small motif to
# identify the most likely 5p end for the 21U-RNA
def predictTou(seq):
    if len(seq)!=1+upLen+downLen:
        raise ValueError("length of sequence was incorrect (" + str(len(seq)) + \
                         " instead of " + str(1+upLen+downLen) + "):\n" + seq + "\n")
    largeScore = largeMotif['apply'](seq[0:largeMotif['len']])
    max_score = None
    max_position = -1
    for n in range(distance['min len'], distance['max len']+1):
        start = largeMotif['len'] + n
        end = start + smallMotif['len']
        new_score = largeScore + distance['apply'](n) + smallMotif['apply'](seq[start:end])
        if max_score==None or new_score>max_score:
            max_score = new_score
            max_position = start + 4  ### this is because the small motif includes four positions
                                      ### upstream of the 21U
    if max_position < 0: raise ValueError("something didn't work; max position is "+str(max_position))
    return max_score,max_position


    
# predicts 21Us along the entire sequence (note: 21U itself must be included in the sequence)
# returns a list of 21U start positions in the sequence (indexed from 0)
def findTous(seq,minScore):
    touStarts = []
    # get all small motif scores ahead of time
    smallScores = map(lambda n: smallMotif['apply'](seq[n:n+smallMotif['len']]), xrange(len(seq) - smallMotif['len']))
    
    # try each big score
    for n in xrange(len(seq) - largeMotif['len'] - smallMotif['len']):
        bigScore = largeMotif['apply'](seq[n:n+largeMotif['len']])
        for d in range(distance['min len'], distance['max len'] + 1):

            # the score must meet the minScore and the 21U must be included in the given sequence
            pos = n + largeMotif['len'] + d
            if pos < len(smallScores) and \
               bigScore + smallScores[pos] + distance['apply'](d) >= minScore and \
               pos + 4 + 21 <= len(seq):
                
                touStarts.append(pos+4)
    return touStarts

                
