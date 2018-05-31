""" 
    This is an implmeentation of a neighbor joining tree based on sequencing
    Right now, it is going to implement based on the rpi lecture 
    cite: http://www.bioinfo.rpi.edu/bystrc/courses/biol4540/lecture.pdf

    Input: none, {soon} it will take in a file or multiple files and do the neighbor joinning tree on it 
    output: tree.txt, {soon} {$input_name}.txt 

    Algorithm input: n x n distance matrix d 
              output: weighted tree t with n leaves fitting d
    
    TODO: figure this part out as well.
    Algorithm: 
    1. make a distance matrix 
    2. 


    Disclaimer: I have not yet ran this program so make sure to debug this once you get the information

"""
import sys, os


# creating the distance-based matrix
# since they're size n x n then we will use a same iteration to figure out the ratio 
# input is the dictionary of sequences 
# output is the ratio of n keys to n keys (as values but in ratio)
# information is based on http://www.bioinfo.rpi.edu/bystrc/courses/biol4540/lecture.pdf
# TODO: this should be right, the answer is on the url link
# TODO: OPTIMIZATION - just compare everything after the key as the ones behind it will be compared already
#                      this should settle the issue if it gets too bad to half of time namely
def creatingRatio(sample): 
    # # honestly being lazy just making this as easy and long as possible LOL 
    # newSample = dict()
    # # initializing the dictionary here 
    # for k, v in sample.iteritems():
    #     newSample[k] = list()
    # # reading each key, then creating the ratio of the letters
    # for k, v in sample.iteritems(): 
    #     # comparing each sample[key] to the others as well
    #     for key in sample: 
    #         if k == key: 
    #             # setting a None for all matching keys
    #             newSample[k].append(None)
    #         else: 
    #             # now comparing the sequences
    #             result = None
    #             for i in range(0,len(v)):
    #                 if v[i] != (newSample[key])[i]:
    #                     result += 1
                
    #             # result = the ratio of non-matching / length of sequence
    #             result = (result / len(v))  
    #             newSample[k].append(result)
    # return newSample
    return 0


# this is the algorithm to work on neighbor joining tree 
# which clustering am I going to do? 
# TODO: make sure you understand what type of clustering algorithm we are settling with 
# here

# I believe we are going to settle with a weighted tree graph and not an unweighed version
def neighborJoining(sample):
    return 0

# step 2 of the algorithm 
# calculate a new distance matrix using each pair 
# M(ij) = d(ij) - [r(i) + r(j)] / (n - 2) 
# ex. M(AB) = d(AB) - [r(A) + r(B)] / (n - 2) = -13 
# or == M(AB) = 5 - [30 + 42] / (6 - 2) == -13
# return a new dictionary with the negative values
def distanceMatrix(targetList): 
    for key, value in targetList.items():
        for i in range(0, len(value)): 
            
# step 1 of the algorithm
# as the name mentions, it is to convert the divergence of each key 
def convertR2Divergence(targetList): 
    # exactly the same order in which the dictionary has the keys as 
    divergence = dict()
    for key, value in targetList.items(): 
        divergence[key] = (sum(value))
    return divergence

# this prints out the current state of the matrices
# prints them out in tab delimited form to make it easier to read
# sorry for the folks using this above 20...
def printMatrices(targetList): 
    n = len(targetList)
    # setting up the columns 
    column = "\t"
    for key in targetList.keys(): 
        column += "{}\t".format(key) 
    print(column)
    currentRow = "" 
    for key, value in targetList.items():
        currentRow = "{}\t".format(key)
        # making this more readable
        tempValue = '\t'.join(str(e) for e in value) 
        currentRow += tempValue
        print(currentRow) 

    # fixing the spacing here
    print("\n")                
def main():

    sample = {'A': "TTGACCAGACCTGTGGTCCG",
              'B': "TTGAACAGACCTGCGGTCGG", 
              'C': "TAGAAAAGACCTGTCGTAGG",
              'D': "GTGCAAAGTCCTGTGTATCG"} 
    
    # assuming this is right
    # cite: http://www.deduveinstitute.be/~opperd/private/neighbor.html
    sampleRatio = {'A': [0, 5, 4, 7, 6, 8],
                   'B': [5, 0, 7, 10, 9, 11],
                   'C': [4, 7, 0, 7, 6, 8], 
                   'D': [7, 10, 7, 0, 5, 9],
                   'E': [6, 9, 6, 5, 0, 8],
                   'F': [8, 11, 8, 9, 8, 0]}
    if sys.stdin: 
        sample = sys.stdin
    
    # printMatrices(sampleRatio)
    neighborJoining(sampleRatio)

if __name__ == "__main__":
    main()


"""
Documentation down here: 
1. http://www.deduveinstitute.be/~opperd/private/neighbor.html for how it works in example 
2. 
"""