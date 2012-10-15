from random import getrandbits

for B in range(1, 6):
    for L in range(5, 21, 5):
        BITS = B * 8
        LINES = L
        print BITS, L
        filename = "generic_lines_" + str(LINES) + "_bits_" + str(BITS) + '.dat'
        f = open(filename, "w")
        f.write(str(LINES) + '\n')

        for j in range(0, L):
            for i in range(0, 6):
                f.write(str(getrandbits(BITS)) + " ")
            f.write('\n')
        
        f.close()


