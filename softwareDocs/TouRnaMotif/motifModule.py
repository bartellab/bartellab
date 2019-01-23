# This module is for building position-specific scoring matrices that
# implement the motif interface.  That interface is:

# Motif Interface
#
# motif = dictionary, with keys representing specfields
#
# motif['name'] = a String of the motif's name
# motif['string'] = a String that describes the motif
# motif['len'] = the int length of the motif, i.e. the number of positions
#     between and including the first position of the motif and the last
#     position of the motif.
# motif['apply'] = a lambda that takes one argument (a string) and returns
#     a number (i.e. a score) indicating how well the input sequence matches
#     the motif.

import math

#############################################################################
######                          SEQUENCE                               ######
#############################################################################

## SEQUENCE is the most basic definition of a sequence motif.  It defines a
## motif as a string of characters.  Matches to the motif are matches to that
## string

## The SEQUENCE constructor takes a string as its argument and returns a motif
## interface which considers the motif to be a perfect match to that sequence.
## @requires sequence is a String

def SEQUENCE(sequence):
    md = dict()
    md['name'] = sequence
    md['string'] = sequence
    md['len'] = len(sequence)
    def make_apply(sequence):
        ## apply returns a 1 if seq is a perfect match to sequence,
        ## otherwise a zero
        ## @requires seq is a String
        def apply(seq):
            if seq == sequence: return 1
            else: return 0
        return apply
    md['apply'] = make_apply(sequence)
    return md

#############################################################################
######                          SEQ_LIST                               ######
#############################################################################

## SEQ_LIST is a definition of a motif where each position of the motif can
## have one in a set of identities.  

## The SEQ_LIST constructor takes a List of Strings as its argument.  Each string
## in the list represents a position in the motif, and each character in the string
## represents an acceptable match for a character in that position.  The special
## character '*' anywhere in the position's string represents that ANY character is
## acceptable in that position.  The constructor returns a motif interface reflecting
## these properties.
## @requires m is a list of strings.  An example of a motif represented thus:
## ['C','T','G','T','T','T','C','A','*','*','*','*','G']

def SEQ_LIST(m):
    ## md is the motif dictionary interface; filled with name, string, len
    md = dict()
    md['name'] = ','.join(m)
    md['string'] = ','.join(m)
    md['len'] = len(m)

    ## make_apply generates apply
    def make_apply(m):
        ## I make a defensive copy of m so that the client can't mutate it to
        ## change apply's behavior.
        my_m = []
        for i in m: my_m.append(i)
        mlen = len(my_m)
        
        ## apply returns a 1 if all of the positions in seq a) match a character from that position
        ## in my_m, or b) the string for that position in my_m contains a '*" character; otherwise
        ## returns 0
        ## @requires seq is a string
        ## @raises ValueError if length of sequence provided doesn't match length of motif
        def apply(seq):
            if len(seq)!=mlen:
                raise ValueError, 'Length of sequence provided ('+str(len(seq))+') is not equal to length of motif ('+str(len_of_motif)+').'
            result = True
            place = 0
            while result and place < mlen:
                check = False
                for c in my_m[place]: check = c=='*' or c==seq[place]
                result = check
                place += 1
            if result: return 1
            else: return 0
        return apply
    md['apply'] = make_apply(m)
    return md



#############################################################################
######                             PSSM                                ######
#############################################################################

# the PSSM constructor will take in a text string which is a table of space- or tab-separated
# columns depicting nucleotide frequencies and construct a PSSM from it.  the format for the
# text string will be like the following example (minus the comment hashes):
# """
# background   .25  .25 .25 .25
# table-name    A   T   G   C
# 0             50 97   200 31
# 5             12  4   350 7
# """

# Requirements of the input string:
#        - lines must be separated by '\n'
#        - first line with text must be background frequencies in float notation
#        - second line with text must have letter keys
#        - lines without text must have no other characters
#        - lines with text must have the same number of entries
# @returns a motif representing a pssm

def PSSM(text):
    ## courtesy fail-fast checks to some input string specs;
    ## these lines make a list from the input string
    lines = filter(lambda a: a!='', text.split('\n'))
    if len(lines)<2: raise SyntaxError, 'Too few lines in data table.'

    ## parse text into a series of lists: name_key is the name, followed by the keys;
    ## background is the word 'background', followed by the background frequencies (floats);
    ## each list in positions is the position #, followed by the counts for that letter in
    ## that position.  Columns all have the same list indexes.  Courtesy exceptions are
    ## thrown if columns are missing.
    background_done = names_done = False
    positions = []
    for l in lines:
        if not(background_done):
            background = ['background']
            background.extend(map(float, filter(lambda a: a!='', l.split())[1:]))
            background_done = True
        elif not(names_done):
            name_key = filter(lambda a: a!='', l.split())
            names_done = True
            if len(name_key)!=len(background): raise SyntaxError, 'Inappropriate number of columns in line: '+str(len(name_key))+' instead of '+str(len(background))+'\n'+l
        else:
            new_line = map(int, filter(lambda a: a!='', l.split()))
            positions.append(new_line)
            if len(name_key)!=len(new_line): raise SyntaxError, 'Inappropriate number of columns in line: '+str(len(new_line))+' instead of '+str(len(name_key))+'\n'+l

    ## score_rep is a dictionary of dictionaries, where the outer keys are position numbers
    ## and the inner keys are letters; values are the scores for each letter at each position
    score_rep = dict()
    motif_start = min(map(lambda a: a[0], positions))
    for p in positions:
        new_pos = dict()
        for n in range(1,len(p)):
            new_pos[name_key[n]] = PSSM_compute_score(p[n], p[1:], background[n], background[1:])
        score_rep[p[0] - motif_start] = new_pos

    ## The string rep of the motif, string_rep, is a tab-delimited print-out of score_rep
    string_rep = '\t'.join(name_key)+'\n'
    n = 0
    while n < max(score_rep.keys())+1:
        if score_rep.has_key(n):
            string_rep += str(n)+'\t'+'\t'.join(map(lambda k: str(round(score_rep[n][k],4)), name_key[1:]))+'\n'
        n+=1
    
    ## the motif interface is made and returned
    motif = dict()
    motif['name'] = name_key[0]
    motif['len'] = max(score_rep.keys())+1
    motif['string'] = string_rep
    motif['apply'] = PSSM_get_apply_from_dict(score_rep)
    return motif


                        
#############################  PSSM Helper Functions  ################################

## This method generates and returns an 'apply' method
def PSSM_get_apply_from_dict(d):
    len_of_motif = max(d.keys())+1
    ## apply returns the sum of the log-odds score for each scored position in the motif
    ## across the supplied seq.
    ## @requires seq is a string
    ## @raises SyntaxError if length of sequence provided doesn't match length of motif
    ## @raises KeyError if the seq contains a letter not provided as a key in the original
    ##      string representation of the motif
    def apply(seq):
        if len(seq)!=len_of_motif:
            raise ValueError, 'Length of sequence provided ('+str(len(seq))+') is not equal to length of motif ('+str(len_of_motif)+').'
        score = 0
        for n in d.keys(): score += d[n][seq[n]]
        return score
    return apply

## This function computes the log-odds score for an letter in a position.  Pseudocounts
## are incorporated into the returned solution according to the method selected.
def PSSM_compute_score(fore, fore_list, back, back_list):
    ## Only one of the two possible methods for calculating fore_freq below should be used;
    ## the other(s) should be commented out:
    
    ### fore_freq = PSSM_pseudocount_laplace(fore, fore_list)
    fore_freq = PSSM_pseudocount_sqrt(fore, fore_list, back, back_list)
    
    ## the back_freq value is constant across all of the fore_freq methods
    back_freq = float(back)/sum(back_list)
    ## the log-odds (base 2) score is returned
    return math.log(fore_freq/back_freq, 2)

## Laplace's single pseudocount method
def PSSM_pseudocount_laplace(fore, fore_list):
    numerator = fore + 1
    denominator = sum(fore_list) + len(fore_list)
    return float(numerator)/denominator

## sqrt(N) pseudeocount method
def PSSM_pseudocount_sqrt(fore, fore_list, back, back_list):
    rooted_counts = math.sqrt(sum(fore_list))
    pseudo_counts = rooted_counts * float(back) / sum(back_list)
    numerator = fore + pseudo_counts
    denominator = sum(fore_list) + rooted_counts
    return float(numerator)/denominator



#############################################################################
######                           TWO_PSSM                              ######
#############################################################################

# TWO_PSSM integrates two PSSMs together, and combines their values with a third
# value, also derived from a log-odds score; that is the distance between the
# two PSSMS.  The starting position of the sequence read is defined by the user
# input.  The shorter of the two PSSMs is slid through the permissible range of
# positions of the input sequence; at each position, its distance score and motif
# score are summed, and the maximum across the set of permitted input values is
# returned as a sum with the score from the other, static PSSM.
#
# the constructor takes in two PSSMs.  these can be generated using the PSSM
# constructor.  they are given in the order in which they are expected to appear
# (5' first). the remaining input is a string reflecting the length distribution.
# its format is similar to that of the PSSM input string, but with integers representing
# distances between the two PSSMs as opposed to letters.  the name in that string
# will be the name given to the motif.  an example is shown below:
# """
# name   3   4    5   6     7   8    10  11
# blah  47  289 5935 4443 2323 1119 334  23
# """
# ther is no background line at the top because the background assumption is an even
# distribution over all available permissable positions.  the string 'blah' can be
# anything but whitespace; it is just a placeholder.  the final parameter for the
# TWO_PSSM constructor is the separator string, just like in the PSSM constructor.
#
# TWO_PSSMs have an additional attribute: two_pssm['positions'](seq) will return
# the starting positions for pssmA and pssmB as a two-item list of integers.  The
# integers will be the indexes of the starting positions for the two pssms in the
# input seq, which is a string.
#
# @param pssmA is an instance of PSSM (the 5' PSSM)
# @param pssmB is an instance of PSSM (the 5' PSSM)
# @param distances is a string with distance frequencies as described above
# @requires pssmA and pssmB are PSSMs
# @requires there is at least one distance specified in 'distances'
# @requires that all number values in 'distances' are integers
# @requires String separator is '\t' or ' '

def TWO_PSSM(pssmA,pssmB,distances):
    ## motif is the motif dictionary
    motif = dict()
    
    ## split distances into a list of lists; list 0 is the top line of the input string, list 1
    ## is the bottom line.  the elements are items from that line as strings
    dist_l1 = filter(lambda a: a!='', distances.split('\n'))
    dist_l2 = map(lambda a: filter(lambda b: b!='', a.split()), dist_l1)

    ## throw exceptions for 'distances' iput violations (no distances; inconsistant numbers of
    ## elements on the two lines of the input string):
    if not(len(dist_l2[0]) > 1):
        raise SyntaxError, 'Distance table must contain at least one entry; yours doesn\'t.\n'+dist_l1[0]
    elif not(len(dist_l2[0])==len(dist_l2[1])):
        raise SyntaxError, 'Inappropriate number of columns in distance table: should be '+str(len(dist_l2[0]))+', but there are '+str(len(dist_l2[1]))+'.\n'+dist_l1[1]

    ## obtain name from distances
    motif['name'] = dist_l2[0][0]

    ## dist_dict's keys will be distances (ints); values will be count #'s (also ints)
    ## dist_background's keys will be distances (ints); values will be naive background frequencies;
    ## 1 / the number of candidate positions.
    dist_dict = dict()
    dist_background = dict()
    for n in range(1,len(dist_l2[0])):
        dist_dict[int(dist_l2[0][n])] = int(dist_l2[1][n])
        dist_background[int(dist_l2[0][n])] = 1.0/(len(dist_l2)-1)
    
    ## get the scores for the distance dictionary, put them into dist_scores.
    ## Note: the PSSM helper function PSSM_compute_score is used here, so the pseudocount scheme
    ## used in the construction of the PSSMs will automatically be applied here as well.
    dist_scores = dict()
    for k in dist_dict.keys():
        dist_scores[k] = PSSM_compute_score(dist_dict[k],dist_dict.values(),dist_background[k],dist_background.values())

    ## length is the sum of the two PSSM's lengths plus the max number of poitions in between
    motif['len'] = pssmA['len'] + pssmB['len'] + max(dist_dict.keys())

    ## get the apply and positions methods using TWO_PSSM_get_apply_positions
    apply_positions = TWO_PSSM_get_apply_positions(pssmA,pssmB,dist_scores)
    motif['apply'] = lambda seq: apply_positions(seq)[0]
    motif['positions'] = lambda seq: apply_positions(seq)[1:]

    ## The string rep of the motif, string_rep, is a tab-delimited print-out each pssm (labelled
    ## as 'pssmA' and 'pssmB') and the distance scoring matrix.
    string_rep = 'pssmA:\n'+pssmA['string']
    string_rep += 'pssmB:\n'+pssmB['string']
    string_rep += 'distances:\n'
    string_rep += '\t'.join(dist_l2[0][1:])+'\n'
    string_rep += '\t'.join(map(lambda k: str(round(dist_scores[int(k)],4)), dist_l2[0][1:]))+'\n'

    ## 'string' is set to the string_rep
    motif['string'] = string_rep

    ## motif can be returned containing everything it needs
    return motif

#############################  TWO_PSSM Helper Functions  ################################

## This method generates and returns an 'apply' method
def TWO_PSSM_get_apply_positions(pssmA,pssmB,dist):
    len_of_motif = pssmA['len'] + pssmB['len'] + max(dist.keys())
    max_distance = max(dist.keys())

    ## apply_positions takes in a sequence and returns a list whose items are 0) the best score for
    ## the two-PSSM matrix, 1) the starting index for pssmA in seq, and 2) the starting index for
    ## pssmB in seq.  that score is determined by taking the shorter of the two pssms and sliding it
    ## through the acceptable range of distances from the larger of the two pssms.  this is a performance
    ## choice; since the slid pssm will need to be score many more times than the static pssm, the smaller
    ## pssm should be selected for that task.  it is not necessarily the case that the larger pssm scores
    ## more positions than the shorter pssm (see PSSM specs), but i anticipate that this will generally be
    ## the case.
    ## @requires seq is a string
    ## @raises SyntaxError if length of sequence provided doesn't match length of motif
    ## @raises KeyError if the seq contains a letter not provided as a key in the original
    ##      string representation of the motif
    
    if pssmA['len'] >= pssmB['len']:
        def apply_positions(seq):
            if len(seq)!=len_of_motif:
                raise ValueError, 'Length of sequence provided ('+str(len(seq))+') is not equal to length of motif ('+str(len_of_motif)+').'
            distK = dist.keys()
            indexA = 0
            for k in distK:
                local_score = dist[k] + pssmB['apply'](seq[pssmA['len']+k:pssmA['len']+k+pssmB['len']])
                if k==distK[0] or local_score > score:
                    score = local_score
                    indexB = pssmA['len']+k
            score += pssmA['apply'](seq[:pssmA['len']])
            return [score,indexA,indexB]
    else:
        def apply_positions(seq):
            if len(seq)!=len_of_motif:
                raise ValueError, 'Length of sequence provided ('+str(len(seq))+') is not equal to length of motif ('+str(len_of_motif)+').'
            distK = dist.keys()
            for k in distK:
                local_score = dist[k] + pssmA['apply'](seq[max_distance-k:max_distance-k+pssmA['len']])
                if n==0 or local_score > score:
                    score = local_score
                    indexA = max_distance-k
            indexB = max_distance + pssmA['len']
            score += pssmB['apply'](seq[indexB:])
            return [score,indexA,indexB]
    return apply_positions







#############################################################################
######                         DISTANCE_PSSM                           ######
#############################################################################

# This is like the distance portion of TWO_PSSM, but alone.
#
# the constructor takes in a string reflecting the length distribution.
# its format is similar to that of the PSSM input string, but with integers representing
# distances between the two PSSMs as opposed to letters.  the name in that string
# will be the name given to the motif.  an example is shown below:
# """
# name   3   4    5   6     7   8    10  11
# blah  47  289 5935 4443 2323 1119 334  23
# """
# ther is no background line at the top because the background assumption is an even
# distribution over all available permissable positions.  the string 'blah' can be
# anything but whitespace; it is just a placeholder.  the final parameter for the
# constructor is the separator string, just like in the PSSM constructor.
#
# the distance pssm does not fully implement the PSSM interface because the length value
# is meaningless here.  instead, there are two keys: 'min len' and 'max len',
# corresponding to the minimum and maximum key values.

def DISTANCE_PSSM(distances):
    
    ## motif is the motif dictionary
    motif = dict()

    ## split distances into a list of lists; list 0 is the top line of the input string, list 1
    ## is the bottom line.  the elements are items from that line as strings
    dist_l1 = filter(lambda a: a!='', distances.split('\n'))
    dist_l2 = map(lambda a: filter(lambda b: b!='', a.split()), dist_l1)

    ## throw exceptions for 'distances' iput violations (no distances; inconsistant numbers of
    ## elements on the two lines of the input string):
    if not(len(dist_l2[0]) > 1):
        raise SyntaxError, 'Distance table must contain at least one entry; yours doesn\'t.\n'+dist_l1[0]
    elif not(len(dist_l2[0])==len(dist_l2[1])):
        raise SyntaxError, 'Inappropriate number of columns in distance table: should be '+str(len(dist_l2[0]))+', but there are '+str(len(dist_l2[1]))+'.\n'+dist_l1[1]

    ## obtain name from distances
    motif['name'] = dist_l2[0][0]

    ## dist_dict's keys will be distances (ints); values will be count #'s (also ints)
    ## dist_background's keys will be distances (ints); values will be naive background frequencies;
    ## 1 / the number of candidate positions.
    dist_dict = dict()
    dist_background = dict()
    for n in range(1,len(dist_l2[0])):
        dist_dict[int(dist_l2[0][n])] = int(dist_l2[1][n])
        dist_background[int(dist_l2[0][n])] = 1.0/(len(dist_l2)-1)
    
    ## get the scores for the distance dictionary, put them into dist_scores.
    ## Note: the PSSM helper function PSSM_compute_score is used here, so the pseudocount scheme
    ## used in the construction of the PSSMs will automatically be applied here as well.
    dist_scores = dict()
    for k in dist_dict.keys():
        dist_scores[k] = PSSM_compute_score(dist_dict[k],dist_dict.values(),dist_background[k],dist_background.values())

    # @requires that n is a key in the distance range (inclusive of the ends)
    motif['apply'] = lambda n: dist_scores[n]

    motif['min len'] = min(dist_scores.keys())
    motif['max len'] = max(dist_scores.keys())

    string = ''
    for n in range(min(dist_scores.keys()), max(dist_scores.keys())+1):
        string += str(n)+'\t'+str(dist_scores[n])+'\n'

    motif['string'] = string

    return motif
