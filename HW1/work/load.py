def loadSeq(fileName):
    #open file
    with open(fileName, 'r') as file:
        #[Question] Do all file only have trash info in first line?####################
        #read first line, and ignore it
        file.readline()
        
        #store info without newline
        DNA_info = file.read().replace('\n', '')
    #close file
    file.close()

    return str(DNA_info)
