# -*- coding: utf-8 -*-

import os, sys, re

normal_nts = ['A', 'C', 'G', 'T']

def ATtrick(_bim, _out):

    with open(_bim, 'r') as f_bim, \
            open(_out+'.ATtrick.bim', 'w') as f_out1, \
            open(_out+'.ATtrick.a1_allele', 'w') as f_out2:

        for line in f_bim:

            l = re.split(r'\s+', line.rstrip('\n'))

            """
            0: chr
            1: Label
            2: GD
            3: BP
            4: a1
            5: a2
            """

            a1 = l[4]
            a2 = l[5]

            if (a1 in normal_nts) and (a2 in normal_nts):
                f_out1.write(line)
                f_out2.write('\t'.join([l[1], l[4]]) + "\n")
            else:
                f_out1.write('\t'.join([l[0], l[1], l[2], l[3], 'T', 'A']) + "\n")
                f_out2.write('\t'.join([l[1], 'T']) + "\n")

    return [_out+'.ATtrick.bim', _out+'.ATtrick.a1_allele']




if __name__ == '__main__':

    [_bim, _out] = sys.argv[1:]
    ATtrick(_bim, _out)