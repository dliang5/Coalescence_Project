""" 
    This is an implmeentation of a neighbor joining tree based on sequencing
    Right now, it is going to implement based on the rpi lecture 
    cite: http://www.bioinfo.rpi.edu/bystrc/courses/biol4540/lecture.pdf

    Input: none, {soon} it will take in a file or multiple files and do the neighbor joinning tree on it 
    output: tree.txt, {soon} {$input_name}.txt 

    Algorithm input: n x n distance matrix d 
              output: weighted tree t with n leaves fitting d
    

    Incoming updates after I finish the example algorithm to this: 
    1. change the key names to numbers in string form 
        1. unless I want to do a permutation of letters which is going to be weird after a while 
        2. Rationale: because there could be over 1k individual genomes (ex. Drosophila mel)
    2. optimize the search in the other functions
    3. create a pattern for the list pattern and create a picture of the tree somehow
    4. make sure the branch length is correct with the example online 

    Concerns and Disclaimer: 
    1. This algorithm is very biased towards whoever is the smallest first
        because of this, it is sub-optimal to the actual optimal one as it does not necessarily produce the correct "lineage" which may
        be of some concern
    2. my branch length may be off so I gotta double check on that later on 
"""
import sys, os, math
from copy import deepcopy

# picking the smallest M(ij) and combine them into U 
# in the current given example, it would be M(AB) or M(DE) in this case 
# I will be combining the earliest ones I can find
# this is when im going to change it to half a graph and not the full one
# as it will mess up the tree
# targetList = matrix after the difference between distances 
# priorList = the state before incoming changes here 
# formedBranches = containing branches used and formed already into new neighbors
def neighborJoining(targetList, priorList, divergences, formedBranches): 
    smallest = 1000000000000
    smallestRowPos = 0
    smallestRowPosLet = ""
    smallestColPosLet = ""
    smallestColPos = 0 
    
    for rowIndex, rowKey in enumerate(targetList): 
        for colIndex, colKey in enumerate(targetList): 
            # leaving them at 0.
            # if rowKey == colKey: 
            #     break 

            # storing the information 
            # i should probably go 
            if smallest > targetList[rowKey][colIndex]: 
                smallest = targetList[rowKey][colIndex] 
                smallestRowPos = rowIndex
                smallestRowPosLet = rowKey
                smallestColPosLet = colKey
                smallestColPos = colIndex
    
    # found the position and their keys to combine
    print("combining {} {} and d({}{}) = {}".format(smallestRowPosLet, smallestColPosLet, 
                                                    smallestRowPosLet, smallestColPosLet,
                                                    smallest))
    
    # getting the 2 branch lengths from smallestColPosLet and smallestRowPos to new branch U
    # S(AU) =d(AB) / 2 + [r(A)-r(B)] / 2(N-2) = 1 
    # S(BU) =d(AB) -S(AU) = 4
    # print("{} {} {}".format( priorList[smallestRowPos][smallestColPos] / 2, divergences[smallestColPosLet] - divergences[smallestRowPos],(2* (len(priorList) - 2) ) ))
    branchLen1 = ((priorList[smallestRowPosLet][smallestColPos] / 2) + (divergences[smallestColPosLet] - divergences[smallestRowPosLet]) ) / (2* (len(priorList) - 2) )
    branchLen1 = math.floor(abs(branchLen1))
    branchLen2 = (priorList[smallestRowPosLet][smallestColPos]) - branchLen1

    # print("this is S({}U) = {}, this is S({}U) = {}".format(smallestRowPosLet, branchLen1, 
    #                                                         smallestColPosLet, branchLen2))


    # deleting all information related to the two branches
    # conditional to delete the furtherest one atm  because it will move one back if i delete the earliest on
    newList = deepcopy(priorList)
    firstDel = 0
    secondDel = 0 
    if smallestColPos > smallestRowPos:
        firstDel = smallestColPos
        secondDel = smallestRowPos
    else:
        firstDel = smallestRowPos 
        secondDel = smallestColPos

    for key in newList:
        del newList[key][firstDel]
        del newList[key][secondDel]

    # now the current dictionary is going to be changed here
    del newList[smallestColPosLet]
    del newList[smallestRowPosLet]

    # initalizing a new dictionary and list.
    newNeighborName = smallestColPosLet + smallestRowPosLet

    # creating a list of where all the branches point to
    # [combinedBranch, self.branch, branchLength]
    formedBranches[smallestColPosLet] = [newNeighborName, smallestColPosLet, branchLen2]
    formedBranches[smallestRowPosLet] = [newNeighborName, smallestRowPosLet,branchLen1]

    # newNeighbor[newNeighborName] = [0] * (len(newList) - 2)
    newList[newNeighborName] = list()
    # (unoptimized) modifying the dictionary with the new branch, newNeighbor
    for rowIndex, rowKey in enumerate(newList):
        entry = 0 
        if rowKey != newNeighborName:
            entry = ( (priorList[rowKey][smallestRowPos] + priorList[rowKey][smallestColPos]) - \
                    priorList[smallestRowPosLet][smallestColPos]) / 2

            # print("{} {} {}".format(priorList[rowKey][smallestRowPos], priorList[rowKey][smallestColPos], 
            #                         priorList[smallestRowPosLet][smallestColPos]))

            # print("d({}U) = d({}{}) + d({}{}) - d({}{}) / 2 = {}".format(rowKey, smallestRowPosLet, rowKey,
            #                                                             smallestColPosLet, rowKey, 
            #                                                             smallestRowPosLet, smallestColPosLet,
            #                                                             entry))
        # remmeber TODO: think about the placement of the new branch because I have no idea where LOL 
        # right now it is selected at the back so it's going to be mess as hell 
        newList[newNeighborName].append(entry)
        if rowKey != newNeighborName:
            newList[rowKey].append(entry)

    # printMatrices(newList)
    return newList

# step 2 of the algorithm 
# calculate a new distance matrix using each pair 
# M(ij) = d(ij) - [r(i) + r(j)] / (n - 2) 
# ex. M(AB) = d(AB) - [r(A) + r(B)] / (n - 2) = -13 
# or == M(AB) = 5 - [30 + 42] / (6 - 2) == -13
# return a new dictionary with the negative values
def distanceMatrix(targetList, divergences): 
    currentList = deepcopy(targetList)
    n = len(targetList) 
    for rowIndex, rowKey in enumerate(targetList): 
        # looping through the list, this is where O(n^3) comes in.
        for colIndex, colKey in enumerate(targetList): 
            if rowIndex == colIndex: 
                targetList[rowKey][colIndex] = 0 
                continue

            # assuming this is after diagonal values
            value = targetList[rowKey][colIndex] - ((divergences[rowKey] + divergences[colKey]) / (n-2))  
            targetList[rowKey][colIndex] = value
    return targetList, currentList
            
# step 1 of the algorithm
# as the name mentions, it is to convert the divergence of each key 
def convertR2Divergence(targetList): 
    # exactly the same order in which the dictionary has the keys as 
    divergence = dict()
    for key, value in targetList.items(): 
        divergence[key] = (sum(value))
    return divergence

# creating the distance-based matrix
# since they're size n x n then we will use a same iteration to figure out the ratio 
# input is the dictionary of sequences 
# output is the ratio of n keys to n keys (as values but in ratio)
# information is based on http://www.bioinfo.rpi.edu/bystrc/courses/biol4540/lecture.pdf
# TODO: this should be right, the answer is on the url link
# TODO: OPTIMIZATION - just compare everything after the key as the ones behind it will be compared already
#                      this should settle the issue if it gets too bad to half of time namely
# this is assuming the user returns a sequence matrix
def creatingRatio(sample): 
    # # honestly being lazy just making this as easy and long as possible LOL 
    newSample = dict()
    # initializing the dictionary here 
    for k, v in sample.items():
        newSample[k] = list()

    for key in sample: 
        for secKey in sample: 
            if key == secKey:
                newSample[key].append(0)
            else: 
                result = 0 
                for i in range(0, len(sample[key])): 
                    if sample[key][i] != sample[secKey][i]:
                        result += 1
                result = result / len(sample[key])
                newSample[key].append(result)

    # printMatrices(newSample)
    return newSample

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

# creating the pattern to the branches 
def branchArt(branches): 
    for key in branches: 
        print(key, branches[key])
    
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

    # printMatrices(sampleRatio)
    creatingRatio(sample)

    formedBranches = dict()
    n = len(sampleRatio)
    for i in range(0, n-2):
        divergences = convertR2Divergence(sampleRatio) 
        nextRatio, sampleRatio = distanceMatrix(sampleRatio, divergences)
        sampleRatio = neighborJoining(nextRatio, sampleRatio, divergences, formedBranches)

    branchArt(formedBranches)
if __name__ == "__main__":
    main()


"""
Documentation down here: 
1. http://www.deduveinstitute.be/~opperd/private/neighbor.html for how it works in example 
2. http://www.srmuniv.ac.in/sites/default/files/files/3(5).pdf
"""