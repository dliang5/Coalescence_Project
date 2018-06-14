""" 
    This will take in a VCF and return annotation of where to find the selected details for the breakpoints of 2Rb
    There is a few sections of how it will produce these details
    1. tandem duplication - <DUP:TANDEM> and their respective ranges
    2. Ref Length and Alt length if there is a huge breakpoint 
    3. Inversion with the ] ] and [ [ style 
        1. inversion breakpoints should have ] ] and [ [ in the same area which is a minor guess


    Input: VCF file 
    Output: parsed reads and information about those reads
"""

import sys, os

def main(): 
    with open(sys.argv[1], 'r') as f: 
        tandem = list() 
        inversion = list() 
        for line in f: 
            content = (line.strip("\n")).split("\t")

            # this is looking for tandem duplication 
            if "DUP:TANDEM" in content[4]:
                tandem.append(line) 
            
            if len(content[3]) < len(content[4]): 
                # check for the inversions side 
                if "]" in content or "]" in content: 
                    contentCheck = content[4].replace("]", " ").replace("[", " ").replace(",", " ").split(" ") 
                    if content[0] not in contentCheck: 
                        continue
                    else: 
                        inversion.append(line) 
            
        

            
if __name__ == "__main__":
    main()
