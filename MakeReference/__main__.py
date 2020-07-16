#-*- coding: utf-8 -*-
# python -m MakeReference

# import os, sys, re
import argparse, textwrap

from MakeReference.MakeReference_v2 import MakeReference_v2



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        MakeReference_v2.py

        

    #################################################################################################
                                     '''),
                                     add_help=False,
                                     prog='MakeReference')

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--variants", help="\nInput variants data file(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("--chped", help="\nHLA Type Data(.chped)\n\n", required=True)
    parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg",
                        default="19")
    parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", required=True)
    parser.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", required=True)

    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')
    parser.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')

    # Beagle4.1
    parser.add_argument("--phasing", help="\nPerform phasing with Beagle4.1.\n\n", action='store_true')
    parser.add_argument("--mem", help="\nJava Heap Memory size for Bealge4.1. (ex. 2g, 500m)\n\n", default="4g")
    parser.add_argument("--nthreads", help="\nThe number of threads to use in Bealge4.1. (ex. 2)\n\n", default=1, type=int)

    ##### <for Test> #####

    # 2019. 01. 10
    # args = parser.parse_args(["-chped", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU_HLA.imgt370.4field.chped",
    #                           "-variants", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU",
    #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_2_b_MarkerGenerator/20190110_bMarkerTest/HAPMAP_CEU_HLA.imgt370.hg18",
    #                           "-hg", "18",
    #                           "-dict-AA", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_AA.hg18.imgt370",
    #                           "-dict-SNPS", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_SNPS.hg18.imgt370"
    #                           ])

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    MakeReference_v2(_CHPED=args.chped, _OUT=args.out, _hg=args.hg, _dictionary_AA=args.dict_AA,
                     _dictionary_SNPS=args.dict_SNPS, _variants=args.variants, _java_heap_mem=args.mem,
                     _p_dependency=args.dependency, f_save_intermediates=args.save_intermediates,
                     f_phasing=args.phasing, _nthreads=args.nthreads)