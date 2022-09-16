#
from gettext import find
from re import L
from load import *
from dna import *
import random

# find ORF and store all of them into a list
def oneFrame(DNA):
    # store all ORFs 
    orfList = []
    # rest part of codons which did not start from it
    rest = DNA

    # check whether ORFs' length is < 3, if yes stop now
    if len(DNA) < 3:
        #print("DNA sequence's length is not enough")
        return "DNA sequence's length is not enough"
        ###Another way
        '''
        print("DNA sequence's length is not enough")
        a = []
        return a
        '''      

    # find ORFs(ignore last words which is not enough number of 3)
    for i in range (len(DNA)//3 + len(DNA)%3):
        #store current ORF
        curList = []
        # gain first three using for coding to AA
        cur = rest[0:3]
        # find ORF
        if cur == 'ATG':
            curList.append('ATG')
            #rest of sequence after ATG
            nextRest = rest[3:]
            #[Question]if length is not the multiple of three(when till end), still count all till end?##########################################s
            for j in range ((len(rest)-3)//3 + (len(rest)-3)%3):
                #nextCur is the first three label of nextRest
                nextCur = nextRest[0:3]
                nextRest = nextRest[3:]
                if (nextCur != 'TAA') and (nextCur != 'TAG') and (nextCur != 'TGA'):
                    curList[0] = curList[0] + nextCur
                else:#stop loop
                    break
            #add it to list
            orfList.append(curList[0])                      
        # keep rest sequence except first three
        rest = rest[3:]
    return orfList

#idea from zhanhao
#find only one longest orf
def oneFrameV2(DNA):
    '''
# check whether ORFs' length is < 3, if yes stop now
    if len(DNA) < 3:
        #print("DNA sequence's length is not enough")
        return "DNA sequence's length is not enough"
    '''
    
        
    orfList = []#store final ans
    startIdx = -1
    stopIdx = 0
    for i in range (len(DNA) // 3):
        #get current three label
        cur = DNA[stopIdx:stopIdx+3]      
        if (cur == 'ATG') and (startIdx == -1):
            startIdx = stopIdx
        elif  (cur == 'TAA' or cur == 'TAG' or cur == 'TGA') and (startIdx != -1):
            orfList.append(DNA[startIdx:stopIdx])
            startIdx = -1
        stopIdx = stopIdx + 3
    #consider condition like this: ATG, ATGC, ATGCC, ATGCCCGGG(only ATG, ATG with no stop signal)
    if startIdx != -1:
        orfList.append(DNA[startIdx:])
    return orfList#list 

def longestORF(DNA):
    #check 
    if len(DNA) < 3:
        print("DNA sequence's length is not enough")
        return ""

    orfLongest = ''#store final return in this funciton
    #return a list which length is 3, index 1 represent DNA[0:], 2 represent DNA[1:]
    for i in range (3):
        if len(DNA[i:]) < 3:#check
            break
        curListLongestOrf = [] #list with only one longest ORF(will only be ['xxx'])
        curList = oneFrameV2(DNA[i:])
        if len(curList) > 1:
            #[Question]: what if length is equal? Now if equal, do nothing###################################
            curListLongestOrf.append(curList[0])
            #find largest orf in list
            for i in range (len(curList)):
                if i > 0:
                    if len(curList[i]) > len(curListLongestOrf[0]):
                        curListLongestOrf[0] = curList[i]
        elif len(curList) == 1:
            curListLongestOrf.append(curList[0])
        if len(curListLongestOrf) > 0:#if there are at least one thing in curListLongestOrf, then do compare
            if len(curListLongestOrf[0]) > len(orfLongest):
                orfLongest = curListLongestOrf[0]
    return orfLongest # string

#find longest in both strands
def longestORFBothStrands(DNA):
    if len(DNA) < 3:
        print ("DNA sequence's length is not enough")
        return ""
    
    orfPosLongest = longestORF(DNA)
    orfNegLongest = longestORF(reverseComplement(DNA))
    if len(orfPosLongest) > len(orfNegLongest):
        return orfPosLongest # str
    else:
        return orfNegLongest # str


def longestORFNoncoding(DNA, numReps):# DNA:STR numReps:int
    longestOrfLen = 0
    for i in range (numReps):
        shuffleList = list(DNA)
        random.shuffle(shuffleList) # shuffle it
        shuffleStr = collapse(shuffleList) # collapse shuffled list to a str
        curOrf = longestORFBothStrands(shuffleStr)
        #print(curOrf)
        #print(len(curOrf))
        if len(curOrf) > longestOrfLen:
            longestOrfLen = len(curOrf)

    return longestOrfLen# return a number which idicate the length of longest orf.

#['a','c','b'] -> 'acb', L is list
def collapse(L):
    return ''.join(L) #str


def findORFs(DNA):
    #check 
    if len(DNA) < 3:
        print("DNA sequence's length is not enough")
        return ""

    orfListPos = []#store final return in this funciton
    #return a list which length is 3, index 1 represent DNA[0:], 2 represent DNA[1:]
    for i in range (3):
        if len(DNA[i:]) < 3:#check
            break
        curList = oneFrameV2(DNA[i:])
        if len(curList) > 1:           
            #find largest orf in list
            for i in range (len(curList)):
                #[Question]: what if length is equal? Now if equal, do nothing###################################
                orfListPos.append(curList[i])
        elif len(curList) == 1:
            orfListPos.append(curList[0])
    return orfListPos # list

def findORFsBothStrands(DNA):
    #check 
    if len(DNA) < 3:
        print("DNA sequence's length is not enough")
        return ""

    orfListBoth = findORFs(DNA)#store positive first
    orfListNeg = findORFs(reverseComplement(DNA))
    if len(orfListNeg) > 0:#if orfListNeg have at least one item
        for i in range (len(orfListNeg)):
            orfListBoth.append(orfListNeg[i])
    return orfListBoth#list


def getCoordinates(orf,DNA):
    startIdx = -1
    endIdx = -1

    startIdx = DNA.find(orf)
    if startIdx == -1:
        startIdx = DNA.find(reverseComplement(orf))

    if startIdx != -1:
        endIdx = startIdx + len(orf)

    orfCoordinates = [startIdx,endIdx]
    return orfCoordinates


def geneFinder(DNA, minLen):
    orfList = findORFsBothStrands(DNA)
    orfListBigger = []
    for i in range (len(orfList)):
        if len(orfList[i]) > minLen:#[Question]longer or equal to[>/>=]?######################################
            orfListBigger.append(orfList[i])
    
    finalOutputList = [] 
    for i in range (len(orfListBigger)):
        curOrfInfo = getCoordinates(orfListBigger[i],DNA)
        curOrfInfo.append(codingStrandToAA(orfListBigger[i]))
        finalOutputList.append(curOrfInfo)
    #[Question]is this engough for sort??
    finalOutputList.sort()
    return finalOutputList #list

#1.改进
#2.把<3情况酌情删除
def printGenes(geneList):
    outputInfo = ''
    #[Question] Is this format good?###########
    for i in range (len(geneList)):
        print("startCoord: "+ str(geneList[i][0]))
        print("endCoord: "+ str(geneList[i][1]))
        print("proteinSequence: "+ str(geneList[i][2]) + "\n")


    f = open("GeneFinder.txt","w+")
    #[Question] Is this format good?###########
    for i in range (len(geneList)):
        f.write("startCoord: "+ str(geneList[i][0]) + "\n")
        f.write("endCoord: "+ str(geneList[i][1]) + "\n")
        f.write("proteinSequence: "+ str(geneList[i][2]) + "\n\n")
    f.close()
    #[Question]Do I need return something?#################################
    return "Done"


#########################delete###############
if __name__ == "__main__":
    #X73525 = loadSeq("X73525.fa")
    #thresholdLen = longestORFNoncoding(X73525,1500)
    #geneList1 = geneFinder(X73525, 500)
    #printGenes(geneList1)
    #print(geneList1)

    '''
    a = geneFinder('ATGCCCTGAATGCCCGGGTGAATGCCCGGGCCC',6)
    print(a)
    '''

    '''
    a =  getCoordinates("GTT", "ACGTTCGA")
    print(a)
    b =  getCoordinates("CGAA", "ACGTTCGA")
    print(b)
    c =  getCoordinates("CCAA", "ACGTTCGA")
    print(c)
    '''
    
    
    '''
    a = findORFsBothStrands('ATGAAACAT')
    print(a)
    b =  findORFsBothStrands('ATG')
    print(b)
    c =  findORFsBothStrands('')
    print(c)
    d =  findORFsBothStrands('ATGAAA')
    print(d)
    '''
    
    '''
    a =  findORFs("ATGGGATGAATTAACCATGCCCTAA")
    print(a)
    b = findORFs("GGAGTAAGGGGG")
    print(b)
    c = findORFs("ATG")
    print(c)
    d = findORFs("ATGC")
    print(d)
    d = findORFs("ATGCCCGGG")
    print(d)
    d = findORFs("A")
    print(d)
    '''
    
    '''
    a =  longestORFNoncoding('CT',5)
    print("longest is: " + str(a))

    X73525 = loadSeq("X73525.fa")
    a =  longestORFNoncoding(X73525,50)
    print("longest is: " + str(a))

    '''

    '''
    a = "abc"
    print(a[0:])
    print(a[1:])
    '''
   
    '''
    a =  longestORFBothStrands('CTATTTCATG')
    print(a)

    b =  longestORFBothStrands('C')
    print(b)

    c = longestORFBothStrands('ATG')
    print(c)
    c = longestORFBothStrands('ATGGGGCCC')
    print(c)
    c = longestORFBothStrands('CCAT')
    print(c)   
    '''

    '''
    a =  longestORF('ATGAAATAG')
    print(a)

    b =  longestORF('CATGAATAGGCCCA')
    print(b)

    c =  longestORF('CTGTAA')
    print(c)

    d = longestORF('ATGCCCTAACATGAAAATGACTTAGG')
    print(d)
    
    e = longestORF('C')
    print(e)
    
    e = longestORF('GCAT')#############################?#############################
    print(e)
    '''
    a = oneFrameV2("CCCATGTTTTGAAAAATGCCCGGGTAAA")
    print(a)

    b = oneFrameV2("ATGCCCGGG")
    print(b)

    c =  oneFrameV2("CCATGTAGAAATGCCC")
    print(c)

    d =  oneFrameV2("ATGCCCATGGGGAAATTTTGACCC")
    print(d)

    e = oneFrameV2("ATGCCCCGGG")
    print(e)
    e = oneFrameV2("ATG")
    print(e)
    e = oneFrameV2("AT")
    print(e)
    
    '''
    a = oneFrameV2("CCCATGTTTTGAAAAATGCCCGGGTAAA")
    print(a)

    b = oneFrameV2("ATGCCCGGG")
    print(b)

    c =  oneFrameV2("CCATGTAGAAATGCCC")
    print(c)

    d =  oneFrameV2("ATGCCCATGGGGAAATTTTGACCC")
    print(d)

    e = oneFrameV2("ATGCCCCGGG")
    print(e)
    '''


    '''
    a = oneFrame("CCCATGTTTTGAAAAATGCCCGGGTAAA")
    print(a)

    b = oneFrame("ATG")
    print(b)

    c =  oneFrame("CCATGTAGAAATGCCC")
    print(c)

    d =  oneFrame("ATGCCCATGGGGAAATTTTGACCC")
    print(d)

    e = oneFrame("ATG")
    print(e)
    '''



