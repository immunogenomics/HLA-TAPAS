#-*- coding: utf-8 -*-
# python -m SNP2HLA

# import os, sys, re
import argparse, textwrap
from SNP2HLA.SNP2HLA import SNP2HLA



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='SNP2HLA',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     add_help=False,
                                     description=textwrap.dedent('''\
    #################################################################################################

        < SNP2HLA.py >

        SNP2HLA: Imputation of HLA amino acids and classical alleles from SNP genotypes

        Author: Sherman Jia (xiaomingjia@gmail.com)
                + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
                + Extensive modifications by Phil Stuart (pstuart@umich.edu) on 4/19/16 to allow use of Beagle 4.1:
                    verified to work with 22Feb16, 15Apr16, and 03May16 versions of Beagle.
                + Small modifications by Yang Luo (yangluo@broadinstitute.org): 09/30/16: verified working with Bealge 4.1 27Jun16.
                + Recoded to Python and updated by Wanson Choi(wansonchoi@snu.ac.kr) : 2019/02/06, 2020/04/30


        DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.        

        INPUTS:
        1. Plink dataset (*.bed/bim/fam)
        2. Reference dataset (*.bgl.phased.vcf.gz(Beagle 4.1))

        DEPENDENCIES: (download and place them in the 'dependency/' folder.)
        1. PLINK (v1.9)  (Will not work with older Plink 1.07)
        2. Beagle (v4.1) (Need to rename java executable as 'beagle.jar')
        3. vcf2gprobs.jar (Beagle utility for generating a Beagle v3 genotypes probability file from a Beagle 4.1 vcf file with GT field data)
        4. [Optional] If genetic_map_file argument is specified, PLINK format genetic map on cM scale 
            (plink.chr6.GRCh36.map, downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)


    #################################################################################################
                                     '''))

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--target", "-t", help="\nThe prefix of Target PLINK SNP data to impute.(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)
    parser.add_argument("--reference", "-ref", help="\nThe prefix of reference panel for imputation.\n\n", required=True)

    parser.add_argument("--tolerated-diff", help="\nTolerated diff (default : 0.15).\n\n", default=0.15)
    parser.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')

    # Beagle(v4).
    parser.add_argument("--mem", help="\nJava Heap Memory size for Bealge4.1. (ex. 2g, 500m)\n\n", default="4g")
    parser.add_argument("--nthreads", help="\nThe number of threads to use in Bealge4.1. (ex. 2)\n\n", default=1, type=int)
    parser.add_argument("--niterations", help="\n'niterations' parameter in Beagle4.1 (default: 5).\n\n", default=5)

    ##### <for Test> #####

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    SNP2HLA(args.target, args.reference, args.out, _mem=args.mem, _tolerated_diff=args.tolerated_diff,
            _p_dependency=args.dependency, _b_nthreads=args.nthreads, _b_niterations=args.niterations)

