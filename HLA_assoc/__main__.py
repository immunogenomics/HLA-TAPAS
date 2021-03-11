# -*- coding: utf-8 -*-
# python -m HLA_assoc

#import os, sys, re
import argparse, textwrap

from HLAassoc.HLAassoc import HLAassoc



if __name__ == "__main__":


    parser = argparse.ArgumentParser(prog='HLAassoc',
                                     add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        HLAassoc

        Supervisor : WansonChoi
        E-mail : wansonchoi@snu.ac.kr

        Perform single marker logistic/linear regression using PLINK.

    #################################################################################################
                                     '''))


    parser._optionals.title = "OPTIONS"
    parser.add_argument("--help", "-h", help="Show this help message and exit.", action='help')
    # parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n")


    # Subparsers.
    subparsers = parser.add_subparsers(help='Methods for Analysis', dest='Main_Menu')


    ### subparser 1 : Logistic Regression
    subp1 = subparsers.add_parser('LOGISTIC',
                                  help='Logistic Regression',
                                  add_help=False,
                                  formatter_class=argparse.RawTextHelpFormatter,
                                  description=textwrap.dedent('''\
    #################################################################################################
    
        Logistic Regression
        
        Supervisor : WansonChoi
        E-mail : wansonchoi@snu.ac.kr
        
        Description here...
    
    #################################################################################################
    '''))

    subp1._optionals.title = "OPTIONS"

    subp1.add_argument("--help", "-h", help="\nShow this help message and exit.\n\n", action='help')

    subp1.add_argument("--vcf", help="\nIMPUTED or PHASED vcf file to perform logistic regression. (*.vcf(.gz))\n\n", required=True)
    subp1.add_argument("--out", "-o", help="\nOutput file name.\n\n", required=True)
    subp1.add_argument("--reference-bim", help="\nReference bim file to decode ATtrick.\n\n")
    subp1.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')

    subp1.add_argument("--covar", help="\nSpecify .covar file (Plink v1.9).\n\n")
    subp1.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use. (Plink v1.9)\n\n")

    subp1.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.9).\n\n")
    subp1.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use. (Plink v1.9)\n\n")

    CondVars = subp1.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify a single variant ID to condition(i.e. To set it as covariate). (Plink v1.9)\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify a tsv file of multiple variant IDs to condition(i.e. To set it as covariates). (Plink v1.9)\n\n")


    # Reverse-map
    subp1.add_argument("--hped", help="\nHLA Type Data not yet processed by \'NomenCleaner\'.\n\n")
    subp1.add_argument("--chped", help="\nHLA Type Data processed by \'NomenCleaner\' (*.chped)\n\n")
    subp1.add_argument("--hat", help="\nHLA Allele Table file(*.hat).\n\n")



    ### subparser 3 : Linear Regression for continuous trait
    subp3 = subparsers.add_parser('LINEAR',
                                  help='Linear Regression',
                                  add_help=False,
                                  formatter_class=argparse.RawTextHelpFormatter,
                                  description=textwrap.dedent('''\
    #################################################################################################
    
        Performing Linear Regression using PLINK
        
        
    #################################################################################################
    '''))

    subp3._optionals.title = "OPTIONS"

    subp3.add_argument("--help", "-h", help="\nShow this help message and exit.\n\n", action='help')

    subp3.add_argument("--vcf", help="\nIMPUTED or PHASED vcf file to perform logistic regression. (*.vcf(.gz))\n\n", required=True)
    subp3.add_argument("--out", "-o", help="\nOutput file name.\n\n", required=True)
    subp3.add_argument("--reference-bim", help="\nReference bim file to decode ATtrick.\n\n")
    subp3.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')

    subp3.add_argument("--covar", help="\nSpecify .covar file (Plink v1.9).\n\n")
    subp3.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use. (Plink v1.9)\n\n")

    subp3.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.9).\n\n")
    subp3.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use. (Plink v1.9)\n\n")

    CondVars = subp3.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify a single variant ID to condition(i.e. To set it as covariate). (Plink v1.9)\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify a tsv file of multiple variant IDs to condition(i.e. To set it as covariates). (Plink v1.9)\n\n")


    # Reverse-map if want to reverse back to a certain resolution, e.g. G-group
    subp3.add_argument("--hped", help="\nHLA Type Data not yet processed by \'NomenCleaner\'.\n\n")
    subp3.add_argument("--chped", help="\nHLA Type Data processed by \'NomenCleaner\' (*.chped)\n\n")
    subp3.add_argument("--hat", help="\nHLA Allele Table file(*.hat).\n\n")



    ### subparser 2 : Omnibus Test
    subp2 = subparsers.add_parser('OMNIBUS',
                                  help='Omnibus Test',
                                  add_help=False,
                                  formatter_class=argparse.RawTextHelpFormatter,
                                  description=textwrap.dedent('''\
    #################################################################################################

        Omnibus Test

        Author : Masahiro Kanai; Yang Luo
        E-mail : mkanai@broadinstitute.org ; yangluo@broadinstitute.org

        Performing Omnibus test on the haplotype-level formed by amino acid positions

    #################################################################################################
    '''))

    subp2._optionals.title = "OPTIONS"

    subp2.add_argument("--help", "-h", help="\nShow this help message and exit.\n\n", action='help')

    subp2.add_argument("--vcf", help="\nIMPUTED or PHASED vcf file to perform Omnibus Test. (*.vcf(.gz))\n\n")
    subp2.add_argument("--file", help="\nThe common file prefix to nominate '--phased(*.bgl.phased)', '--fam(*.fam)', "
                                      "'--bim(*.bim)', '--pheno(*.pheno)', '--sex(*.sex)', '--pcs(*.pcs)' arguments "
                                      "all at once.\n\n")
    subp2.add_argument("--pop", help="\nFile name for Population information of samples.\n\n")
    subp2.add_argument("--phased", help="\nFile name for BEAGLE phased information.\n\n")
    subp2.add_argument("--fam", help="\nFile name for PLINK fam file. (*.fam)\n\n")
    subp2.add_argument("--bim", help="\nFile name for PLINK bim file. (*.bim)\n\n")
    subp2.add_argument("--pheno", help="\nFile name for Phenotype information of samples. (Case:1 / Control:0).\n\n")
    subp2.add_argument("--sex", help="\nFile name for Sex information of samples.\n\n")
    subp2.add_argument("--pcs", help="\nFile name for Principal Components information.\n\n")
    subp2.add_argument("--maf-threshold", help="\nThreshold value for Minor Allele Frequency.\n\n", default=0.005)
    subp2.add_argument("--aa-only", help="\nRun association test only for AA changes.\n\n", action='store_true')
    subp2.add_argument("--nthreads", help="\nThe number of threads to be used in the Omnibus Test.\n\n", type=int, default=1)
    subp2.add_argument("--out", "-o", help="\nOutput file name.\n\n", required=True)

    subp2.add_argument("--remove-samples-by-haplo", action='store_true',
                       help="\nRemove samples based on haplotypes constructed using a subset of AAs.\n\n")
    subp2.add_argument("--remove-samples-aa-pattern", type=str, help="\nPrefix of markers to exclude. (ex. AA_B)\n\n")
    subp2.add_argument("--min-haplo-count", type=int, default=10, help="\nThe minimum number of haplotype count.\n\n")

    subp2.add_argument("--condition", type=str, help="\nFile name for conditioning.\n\n")
    subp2.add_argument("--condition-gene", type=str, help="\nGene name for conditioning.\n\n")
    subp2.add_argument("--exclude-composites", action='store_true',
                       help="\nExclude composite amino acids from omnibus test.\n\n")
    subp2.add_argument("--output-composites", action='store_true',
                       help="\nOutput composite amino acids in .haplo.txt.\n\n")

    subp2.add_argument("--exhaustive", action='store_true',
                       help="\nRun exhaustive testing of AA combinations in a single HLA gene.\n\n")
    subp2.add_argument("--exhaustive-aa-pos", type=int,
                       help="\nSpecify a AA pos of the first element(s) of combinations.\n\n")
    subp2.add_argument("--exhaustive-min-aa", type=int, default=2,
                       help="\nMinimum number of AA positions to form a combination.\n\n")
    subp2.add_argument("--exhaustive-max-aa", type=int, default=2,
                       help="\nMaximum number of AA positions to form a combination.\n\n")
    subp2.add_argument("--exhaustive-no-filter", action='store_true',
                       help="\nDon't filter non-increasing AA combinations.\n\n")


    subp2.add_argument("--mem", help="\nJava Heap Memory size for vcf2beagle.jar. (ex. 2g, 500m)\n\n", default="2g")
    subp2.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')




    ##### <for Test> #####

    # args = parser.parse_args("--target /home/wansonchoi/Projects/yang/HLA-TAPAS/tests/PLINK_vcf2plink/IMPUTED.1958BC.bgl.phased.vcf.gz "
    #                           "--out /home/wansonchoi/Projects/yang/HLA-TAPAS/tests/PLINK_vcf2plink/IMPUTED.1958BC.HLA_assoc "
    #                           "--hped /home/wansonchoi/Projects/yang/HLA-TAPAS/tests/PLINK_vcf2plink/IMPUTED.1958BC.imgt3320.Ggroup.chped "
    #                           "--chped /home/wansonchoi/Projects/yang/HLA-TAPAS/tests/PLINK_vcf2plink/IMPUTED.1958BC.imgt3320.Ggroup_to_4field.chped "
    #                           "--pheno /home/wansonchoi/Projects/yang/HLA-TAPAS/tests/PLINK_vcf2plink/IMPUTED.1958BC.phe "
    #                           "--logistic "
    #                          "--dependency ../dependency "
    #                          "--reference-bim /home/wansonchoi/Projects/yang/HLA-TAPAS/tests/MyRef2/HAPMAP_CEU.REF.bglv4.bim".split(' '))

    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    # Main implementation

    if args.Main_Menu == 'LOGISTIC':
        HLAassoc(args.Main_Menu, args.out, args.dependency,
                  _vcf=args.vcf, _reference_bim=args.reference_bim,
                  _covar=args.covar, _covar_name=args.covar_name, _pheno=args.pheno, _pheno_name=args.pheno_name,
                  _condition=args.condition, _condition_list=args.condition_list,
                  _hped=args.hped, _chped=args.chped, _hat=args.hat)

    if args.Main_Menu == 'LINEAR':
        HLAassoc(args.Main_Menu, args.out, args.dependency,
                  _vcf=args.vcf, _reference_bim=args.reference_bim,
                  _covar=args.covar, _covar_name=args.covar_name, _pheno=args.pheno, _pheno_name=args.pheno_name,
                  _condition=args.condition, _condition_list=args.condition_list,
                  _hped=args.hped, _chped=args.chped, _hat=args.hat)

    elif args.Main_Menu == 'OMNIBUS':
        HLAassoc(args.Main_Menu, args.out, args.dependency,
                  _vcf=args.vcf, _file=args.file, _pop=args.pop, _phased=args.phased, _fam=args.fam, _bim=args.bim,
                  _pheno=args.pheno, _sex=args.sex, _pcs=args.pcs, _maf_threshold=args.maf_threshold,
                  f_aa_only=args.aa_only, _nthreads=args.nthreads, f_remove_samples_by_haplo=args.remove_samples_by_haplo,
                  f_remove_samples_aa_pattern=args.remove_samples_aa_pattern, _min_haplo_count=args.min_haplo_count,
                  _condition=args.condition, _condition_gene=args.condition_gene,
                  f_exclude_composites=args.exclude_composites,
                  f_output_composites=args.output_composites, f_exhaustive=args.exhaustive,
                  _exhaustive_aa_pos=args.exhaustive_aa_pos, _exhaustive_min_aa=args.exhaustive_min_aa,
                  _exhaustive_max_aa=args.exhaustive_max_aa, f_exhaustive_no_filter=args.exhaustive_no_filter,
                  _java_heap_mem=args.mem)

    else:
        print("Error")

