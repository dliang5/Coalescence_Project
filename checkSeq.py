"""
    This is to read the selected output files and go wild
    What it does is search for the repeated sequences within the haplotypes or VCF by their alternative and reference sequences 
    around the breakpoint. If there is anything, it would print out the reads around it and tell me exactly if it was the repeat 
    from Lobo 2009 or from Russ.

    doesn't seem to work
"""


import os, sys 

literatureRepeat = "ACTTTTGCGATTGTCGCAAAAACTTCTGCGA"
russRepeat = "TAAGCTTTTGCTTTATAAACAAATCCTACGACCTATATTATTATATGCATACCCAGTTTGGGAAAAAAGTGCCTCTACGCACTTAAAGAAAATCCAAGTGTTCCAAAACAAATTTCTTAAACGAATCCTTAATTTACCACTTTGGCATAGTACTAGTGACATTCATCAAAGAACAAATATGATAAAAATAATACAATATGCTGATCTACTATTAGAAAAGTATAAAACAAATTGTAGGTCCAGCCCGTGGCGATTTGTAAATAGGTTATATTAAAGTAAATAGACTAGAT"

with open(sys.argv[1], 'r') as f: 
    alternative = "" # has to match either literature or russ to work 
    alternative1 = ""
    reference = "" # has to match literature i think? 
    for line in f: 
        content = (line.strip("\n")).split("\t")
        flag1 = True
        flag2 = True

        # skip the lines we don't care about 
        if len(content) < 5: 
            continue 
        
        # checks if they match or not.
        if len(alternative) == len(literatureRepeat): 
            for i in range(0, len(alternative)): 
                if alternative[i] != literatureRepeat[i]:
                    flag1 = False
            if flag1 == True: 
                print(line) 
        if len(alternative1) == len(russRepeat): 
            for i in range(0, len(alternative1)): 
                if alternative1[i] == russRepeat[i]: 
                    flag2 = False
            if flag2 == True
                print(line) 
                
        # this is to check the sequences now and append them and remove the first element 
        if len(alternative) == len(literatureRepeat): 
            temp = alternative[1:] + content[4] 
            alternative = temp 
        else: 
            alternative += content[4] 
        if len(alternative1) == len(russRepeat): 
            temp = alternative1[1:] + content[4] 
            alternative1 = temp 
        else: 
            alternative1 += content[4]
        print(alternative, alternative1)
