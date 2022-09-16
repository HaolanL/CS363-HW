import dna
import load

if __name__ == "__main__":
    #function test
    '''
    # test load.py
    DNAInfomation = load.loadSeq('U81861.fa')
    print("HERE IS:" + DNAInfomation)

    # test dna.py
    rcDNA = dna.reverseComplement(DNAInfomation)
    print("The RC is:" + rcDNA)

    #test for codingStrandToAA(DNA) in dna.py
    a = "AGTCCCGGGTTT"
    aOut = dna.codingStrandToAA(a)

    c = "ATGCAACAGCTCT"
    cOut = dna.codingStrandToAA(c)
    '''

    #test for codingStrandToAA(DNA) in dna.py
    a = "AGTCCCGGGTTT"
    aOut = dna.codingStrandToAA(a)
    print(aOut)

    b = "ATGCAACAGCTC"
    bOut = dna.codingStrandToAA(b)
    print(bOut)

    c = "ATGCAACAGCTCT"
    cOut = dna.codingStrandToAA(c)
    print(cOut)


    # other test
    '''
    count = 0
    for i in range (2):
        count += 1
    print(count)

    a = "abc123"
    b = '|'
    c = 'B'
    print(a)
    print(a+b)
    print(a+c)

    # test
    a = "123"
    print(len(a))
    '''

    
    

    