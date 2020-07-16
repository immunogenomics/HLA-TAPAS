#-*- coding: utf-8 -*-
# python -m NomenCleaner


# import os, sys, re
import argparse, textwrap
from NomenCleaner.NomenCleaner import HATK_NomenCleaner




if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='NomenCleaner',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     add_help=False,
                                     description=textwrap.dedent('''\
    #########################################################################################

        < NomenCleaner.py >
        
        Author : WansonChoi
        E-mail : wansonchoi@snu.ac.kr

        - Transforms *.hped file to *.chped file.
        - *.hat file must be given as input file.




    #########################################################################################
                                     '''))

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    # Input (1) : *.hped file
    parser.add_argument("--hped", help="\nHLA Type Data with raw HLA allele(ex. 0101).\n\n", dest="hped", required=True)

    # Input (2) : *.hat file
    parser.add_argument("--hat", help="\nHLA Allele Table file(*.hat).\n\n",
                        default='NomenCleaner/data/HLA_ALLELE_TABLE.imgt3320.hat')
    # parser.add_argument("--imgt", help="\nSpecifying the IMGT-HLA version.\n\n", required=True)

    # Ouptut Prefix
    parser.add_argument("--out", "-o", help="\nOutput file prefix.\n\n", required=True)
    parser.add_argument("--leave-NotFound", help="\nLeaving HLA alleles which can't be found in given *.hat file(Novel or Erroneous allele) intact.\n\n", action='store_true')

    # Output format
    output_digit_selection = parser.add_mutually_exclusive_group()
    output_digit_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n",
                                        action="store_true", dest="oneF")
    output_digit_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n",
                                        action="store_true", dest="twoF")
    output_digit_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n",
                                        action="store_true", dest="threeF")
    output_digit_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields\n\n",
                                        action="store_true", dest="fourF")
    output_digit_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n",
                                        action="store_true")
    output_digit_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n",
                                        action="store_true")

    # Flag to remove HLA gene caption.
    parser.add_argument("--NoCaption", help="\nMake converted HLA alleles NOT have HLA gene prefix(ex. \"A*\").\n\n",
                        action='store_true')

    ##### <for Test> #####


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    HATK_NomenCleaner(args.hped, args.hat, args.out,
                      __f_NoCaption=args.NoCaption, __leave_NotFound=args.leave_NotFound,
                      __oneF=args.oneF, __twoF=args.twoF, __threeF=args.threeF, __fourF=args.fourF,
                      __Ggroup=args.Ggroup, __Pgroup=args.Pgroup)
