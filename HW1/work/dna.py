import aminoAcids

# Takes a DNA string as input (in 5' to 3'order) and returns the sequence of 
# its complementary DNA strand, also in 5' to 3' order.
def reverseComplement(DNA):
    # dnaRC to store reverseComplement sequence of DNA sequence
    dnaRC = ''
    for each in reversed(DNA):
        if each == 'A':
            dnaRC = dnaRC + 'T'
        elif each == 'T':
            dnaRC = dnaRC + 'A'
        elif each == 'C':
            dnaRC = dnaRC + 'G'
        else:
            dnaRC = dnaRC + 'C'
        
    return str(dnaRC)


# takes a sequence of DNA nucleotides from the coding strand 
# and returns the corresponding amino acids as a string.
def codingStrandToAA(DNA):
    # Check whether it is divisible by 3
    if len(DNA) % 3 != 0:
        ##[Question] do we need '???' for ''?##################################
        print("'DNA sequence not divisible by 3!'")
        #[Question]return none mean this or just nothing?############################################
        return None
    
    # AA to store AA info coding from DNA sequence
    AA = ''
    # rest is current DNA sequence which still not convert to AA form
    rest = DNA
    for i in range (len(DNA)/3):
        # gain first three using for coding to AA
        cur = rest[0:3]
        # keep rest sequence except first three
        rest = rest[3:]

        # convert cur to AA
        for j in range (len(aminoAcids.codons)):
            for k in range (len(aminoAcids.codons[j])):
                if cur == aminoAcids.codons[j][k]:
                    AA = AA + aminoAcids.aa[j]

    #[Question]need this? eg: AGT -> 'AGT'################################################
    ##AA = "'" + AA + "'"
    return AA