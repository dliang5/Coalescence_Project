# this one parse out the region required
import sys, os

def main():
    with open(sys.argv[1], 'r') as f:
        tandem = list()
        inversion = list()
        for line in f:
            content = (line.strip('\n')).split('\t')

            if len(content) < 5:
                continue

            flag = True
            for i in content[1]:
                if not i.isdigit():
                    flag = False
                    break
            if flag == False:
                continue

            if int(content[1]) > (19023283-3400) and int(content[1]) < (19023283+3400):
                print(content)
            if int(content[1]) > (31508351-3400) and int(content[1]) < (31508351+3400):
                print(content)
            if int(content[1]) > (26758676-3400) and int(content[1]) < (26758676+3400):
                print(content)

if __name__ == "__main__":
    main()


