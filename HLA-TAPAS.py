# -*- coding: utf-8 -*-

import os, sys, re
from os.path import exists, join, dirname, basename
import argparse, textwrap

from NomenCleaner.NomenCleaner import HATK_NomenCleaner
from MakeReference.MakeReference_v2 import MakeReference_v2
from SNP2HLA.SNP2HLA import SNP2HLA
from HLA_assoc.HLA_assoc import HLA_assoc
from Manhattan.manhattan import Manhattan

# Patterns to use
p_Ggroup_HLA_allele = re.compile(r'^(\w+\*)?\d{2,3}(:\d{2,3})+G$')


std_MAIN_PROCESS_NAME = "[%s]: " % (basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "[%s::ERROR]: " % (basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "[%s::WARNING]: " % (basename(__file__))



def HLA_TAPAS(_target, _reference, _hg, _out, _hped, _chped=None,
              f_save_intermediates=False, _dependency='./dependency',
              b_mem = '2g', b_nthreads=1, b_niterations=5,
              _hat='NomenCleaner/data/HLA_ALLELE_TABLE.imgt3320.hat',
              f_oneF=False, f_twoF=False, f_threeF=False, f_fourF=False, f_Ggroup=False, f_Pgroup=False, f_NoCaption=False,
              _dict_AA='MakeReference/data/HLA_DICTIONARY_AA.hg19.imgt3320',
              _dict_SNPS='MakeReference/data/HLA_DICTIONARY_SNPS.hg19.imgt3320', _tolerated_diff=0.15,
              _pheno=None, _pheno_name=None, _condition=None, _condition_list=None,
              _reference_bim=None, _covar=None, _covar_name=None,
              _pop=None, _sex=None, _pcs=None, _maf_threshold=None, f_aa_only=False,
              f_remove_samples_by_haplo=False, f_remove_samples_aa_pattern=False, _min_haplo_count=10, _condition_gene=None,
              f_exclude_composites=False, f_output_composites=False, f_exhaustive=False,
              _exhaustive_aa_pos=None, _exhaustive_min_aa=2, _exhaustive_max_aa=2, f_exhaustive_no_filter=False
              ):

    """
    Main execution part of HLA-TAPAS pipeline.

    [1] NomenCleaner
    [2] MakeReference_v2
    [3] SNP2HLA
    [4] HLA_assoc

    """

    print("")

    ### Exception Handling

    # './dependency' folder
    if not exists(_dependency):
        print(std_ERROR_MAIN_PROCESS_NAME + "Folder for dependent software('{}') can't be found. Please check '--dependency' argument again.".format(_dependency))
        sys.exit()


    # Target data check.
    if not exists(_target+'.bed'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Target SNP data can't be found('{}'). Please check '--target/-t' argument again.".format(_target + '.bed'))
        sys.exit()
    if not exists(_target+'.bim'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Target SNP data can't be found('{}'). Please check '--target/-t' argument again.".format(_target + '.bim'))
        sys.exit()
    if not exists(_target+'.fam'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Target SNP data can't be found('{}'). Please check '--target/-t' argument again.".format(_target + '.fam'))
        sys.exit()


    # Reference data check.
    if not exists(_reference+'.bed'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.bed'))
        sys.exit()
    if not exists(_reference+'.bim'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.bim'))
        sys.exit()
    if not exists(_reference+'.fam'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.fam'))
        sys.exit()
    if not exists(_hped):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA type data for generating Reference panel can't be found('{}'). Please check '--hped' argument again.".format(_hped))
        sys.exit()


    # Dictionary
    if _hg == '18':
        _dict_AA = 'MakeReference/data/hg18/HLA_DICTIONARY_AA.hg18.imgt3320'
        _dict_SNPS = 'MakeReference/data/hg18/HLA_DICTIONARY_SNPS.hg18.imgt3320'
    elif _hg == '19':
        _dict_AA = 'MakeReference/data/hg19/HLA_DICTIONARY_AA.hg19.imgt3320'
        _dict_SNPS = 'MakeReference/data/hg19/HLA_DICTIONARY_SNPS.hg19.imgt3320'
    elif _hg == '38':
        _dict_AA = 'MakeReference/data/hg38/HLA_DICTIONARY_AA.hg38.imgt3320'
        _dict_SNPS = 'MakeReference/data/hg38/HLA_DICTIONARY_SNPS.hg38.imgt3320'
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Bizarre Human Genome version is given. Please check '--hg' argument again.")
        sys.exit()


    # Intermediate path check.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(dirname(_out)):
        out_dirname = dirname(_out)
        os.makedirs(out_dirname, exist_ok=True)
    else:
        out_dirname = ''




    ### [1] NomenCleaner

    # HPED G-group allele check
    if CheckGgroupHPED(_hped) == -1:
        sys.exit()

    print(std_MAIN_PROCESS_NAME + "Generating 6-digit(maximum 4field) CHPED with given HPED('{}')." \
          .format(_hped))

    CHPED = HATK_NomenCleaner(_hped, _hat, _out+'.4field', __f_NoCaption=False, __leave_NotFound=True,
                              __oneF=False, __twoF=False, __threeF=False, __fourF=True, __Ggroup=False, __Pgroup=False)
    # '--leave-NotFound' as True by default as Yang had requested.

    CHPED = CHPED.chped
    print(std_MAIN_PROCESS_NAME + "Generated CHPED: '{}'.\n\n".format(CHPED))



    ### [2] MakeReference_v2
    print(std_MAIN_PROCESS_NAME + "Generating Reference panel.")

    REFERENCE = MakeReference_v2(CHPED, _out + '.REF.bglv4', _hg, _dict_AA, _dict_SNPS, _reference,
                                 _java_heap_mem=b_mem, _p_dependency=_dependency,
                                 f_save_intermediates=f_save_intermediates, f_phasing=True, _nthreads=b_nthreads)
    print(std_MAIN_PROCESS_NAME + "Generated Reference panel : '{}'.\n\n".format(REFERENCE))



    ### [3] SNP2HLA
    print(std_MAIN_PROCESS_NAME + "Performing SNP2HLA imputation.")

    IMPUTED = SNP2HLA(_target, REFERENCE, _out+'.IMPUTED',
                      _tolerated_diff=_tolerated_diff, _p_dependency=_dependency,
                      _mem=b_mem, _b_nthreads=b_nthreads, _b_niterations=b_niterations)
    print(std_MAIN_PROCESS_NAME + "Imputed result : '{}'\n\n".format(IMPUTED))



    ### [4] HLA_assoc

    # LOGISTIC REGRESSION.
    ASSOC_STUDY = HLA_assoc('LOGISTIC', _out + '.IMPUTED', _dependency,
                            _vcf=IMPUTED, _reference_bim=_reference_bim, _covar=_covar, _covar_name=_covar_name,
                            _pheno=_pheno, _pheno_name=_pheno_name, _condition=_condition, _condition_list=_condition_list,
                            _hped=_hped, _chped=CHPED)

    LOGISTIC_RESULT = ASSOC_STUDY.assoc_result
    print(std_MAIN_PROCESS_NAME + "Output Logistic Regression result : '{}'.\n\n".format(LOGISTIC_RESULT))

    # OMNIBUS TEST.
    OMNIBUS = HLA_assoc('OMNIBUS', _out + '.IMPUTED.OMNIBUS', _dependency,
                  _vcf=IMPUTED, _file=None, _pop=_pop, _phased=None, _fam=_target+'.fam', _bim=REFERENCE+'.bim',
                  _pheno=_pheno, _sex=_sex, _pcs=_pcs, _maf_threshold=_maf_threshold,
                  f_aa_only=f_aa_only, _nthreads=1, f_remove_samples_by_haplo=f_remove_samples_by_haplo,
                  f_remove_samples_aa_pattern=f_remove_samples_aa_pattern, _min_haplo_count=_min_haplo_count,
                  _condition=_condition, _condition_gene=_condition_gene,
                  f_exclude_composites=f_exclude_composites,
                  f_output_composites=f_output_composites, f_exhaustive=f_exhaustive,
                  _exhaustive_aa_pos=_exhaustive_aa_pos, _exhaustive_min_aa=_exhaustive_min_aa,
                  _exhaustive_max_aa=_exhaustive_max_aa, f_exhaustive_no_filter=f_exhaustive_no_filter,
                  _java_heap_mem=b_mem)

    OMNIBUS_RESULT = OMNIBUS.omnibus_result
    print(std_MAIN_PROCESS_NAME + "Output Omnibus Test result : '{}'.\n\n".format(OMNIBUS_RESULT))


    ### [5] Manhattan Plot
    print(std_MAIN_PROCESS_NAME + "Plotting Manhattan Plot.")

    MANHATTAN_PLOT = Manhattan([LOGISTIC_RESULT], LOGISTIC_RESULT+'.manhattan', _hg)
    print(std_MAIN_PROCESS_NAME + "Manhattan plot result : '{}'.\n\n".format(MANHATTAN_PLOT))



    print(std_MAIN_PROCESS_NAME + "HLA-TAPAS done.")

    return 0



def CheckGgroupHPED(_hped_Ggroup):

    with open(_hped_Ggroup, 'r') as f_hped_Ggroup:

        flag_stop = False
        count = 0

        for line in f_hped_Ggroup:

            l = re.split(r'\s+', line.rstrip('\n'))

            """
            0: FID,
            1: IID,
            2: PID,
            3: MID,
            4: Sex,
            5: Phe,
            6,7 : HLA-A,
            8,9 : HLA-B,
            ...
            20,21 : HLA-DRB1

            Total 22 columns
            """

            for i in range(6, 22, 2):

                # allele1
                if not (l[i] == '0' or bool(p_Ggroup_HLA_allele.match(l[i]))):
                    print(std_ERROR_MAIN_PROCESS_NAME + "Allele '{}' is not in G-group HLA allele name form." \
                          .format(l[i]))
                    flag_stop = True

                # allele2
                if not (l[i + 1] == '0' or bool(p_Ggroup_HLA_allele.match(l[i + 1]))):
                    print(std_ERROR_MAIN_PROCESS_NAME + "Allele '{}' is not in G-group HLA allele name form." \
                          .format(l[i + 1]))
                    flag_stop = True

                if flag_stop:
                    break

            if flag_stop:
                print("Please check again the row {} of HPED file('{}').".format(count, _hped_Ggroup))
                print("(row {}: {})".format(count, l))
                return -1

            count += 1
            # if count > 5: break

    return 0



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='HLA-TAPAS',
                                     add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        HLA-TAPAS

        Author : Yang Luo
        E-mail : yangluo@broadinstitute.org
        
        Description here...

    #################################################################################################
                                     '''))

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    ### General
    parser.add_argument("--target", "-t", help="\nThe prefix of Target PLINK SNP data to impute.(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("--reference", "-ref", help="\nThe prefix of reference panel for imputation.\n\n", required=True)

    parser.add_argument("--out", "-o", help="\nOutput file prefix.\n\n", required=True)
    parser.add_argument("--hg", help="\nHuman Genome version of Reference panel to generate.\n\n", choices=["18", "19", "38"], metavar="hg", required=True)

    parser.add_argument("--hped-Ggroup", help="\nHLA Type Data with G-group HLA allele label.\n\n", dest="hped", required=True)
    # parser.add_argument("--chped", help="\nHLA Type Data(.chped)\n\n")


    parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')
    parser.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')


    # Beagle4.1
    parser.add_argument("--mem", help="\nJava Heap Memory size for Bealge4.1. (ex. 2g, 500m)\n\n", default="4g")
    parser.add_argument("--nthreads", help="\nThe number of threads to use in Bealge4.1. (ex. 2)\n\n", default=1, type=int)
    parser.add_argument("--niterations", help="\n'niterations' parameter in Beagle4.1 (default: 5).\n\n", default=5)



    ### [1] NomenCleaner
    # parser.add_argument("--hat", help="\nHLA Allele Table file(*.hat).\n\n",
    #                     default='NomenCleaner/data/HLA_ALLELE_TABLE.imgt3320.hat')

    # parser.add_argument("--leave-NotFound", help="\nLeaving HLA alleles which can't be found in given *.hat file(Novel or Erroneous allele) intact.\n\n", action='store_true') # True by default in 'HLA-TAPAS'.

    # Output format
    # output_digit_selection = parser.add_mutually_exclusive_group()
    # output_digit_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n",
    #                                     action="store_true", dest="oneF")
    # output_digit_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n",
    #                                     action="store_true", dest="twoF")
    # output_digit_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n",
    #                                     action="store_true", dest="threeF")
    # output_digit_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields\n\n",
    #                                     action="store_true", dest="fourF")
    # output_digit_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n",
    #                                     action="store_true")
    # output_digit_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n",
    #                                     action="store_true")

    # Flag to remove HLA gene caption.
    # parser.add_argument("--NoCaption", help="\nMake converted HLA alleles NOT have HLA gene prefix(ex. \"A*\").\n\n",
    #                     action='store_true')



    ### [2] MakeReference_v2
    # parser.add_argument("--phasing", help="\nPerform phasing with Beagle4.1.\n\n", action='store_true') # True by default in HLA-TAPAS.
    # parser.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n",
    #                     default='MakeReference/data/HLA_DICTIONARY_AA.hg19.imgt3320')
    # parser.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n",
    #                     default='MakeReference/data/HLA_DICTIONARY_SNPS.hg19.imgt3320')



    ### [3] SNP2HLA
    parser.add_argument("--tolerated-diff", help="\nTolerated diff (default : 0.15).\n\n", default=0.15, type=float)



    ### [4] HLA_assoc

    # Common
    parser.add_argument("--pheno", help="\nSpecify phenotype information file (Plink v1.9).\n\n")
    parser.add_argument("--pheno-name", help="\nSpecify the column name in phenotype file which you will use. (Plink v1.9)\n\n")

    CondVars = parser.add_mutually_exclusive_group()
    CondVars.add_argument("--condition", help="\nSpecify a single variant ID to condition(i.e. To set it as covariate). (Plink v1.9)\n\n")
    CondVars.add_argument("--condition-list", help="\nSpecify a tsv file of multiple variant IDs to condition(i.e. To set it as covariates). (Plink v1.9)\n\n")


    # Logistic Regression
    parser.add_argument("--reference-bim", help="\nReference bim file to decode ATtrick.\n\n")

    parser.add_argument("--covar", help="\nSpecify .covar file (Plink v1.9).\n\n")
    parser.add_argument("--covar-name", help="\nSpecify the column name(s) in .covar file which you will use. (Plink v1.9)\n\n")


    # Omnibus Test
    parser.add_argument("--pop", help="\nFile name for Population information of samples.\n\n")
    parser.add_argument("--sex", help="\nFile name for Sex information of samples.\n\n")
    parser.add_argument("--pcs", help="\nFile name for Principal Components information.\n\n")

    parser.add_argument("--maf-threshold", help="\nThreshold value for Minor Allele Frequency.\n\n", default=0.005)
    parser.add_argument("--aa-only", help="\nRun association test only for AA changes.\n\n", action='store_true')

    parser.add_argument("--remove-samples-by-haplo", action='store_true',
                        help="\nRemove samples based on haplotypes constructed using a subset of AAs.\n\n")
    parser.add_argument("--remove-samples-aa-pattern", type=str, help="\nPrefix of markers to exclude. (ex. AA_B)\n\n")
    parser.add_argument("--min-haplo-count", type=int, default=10, help="\nThe minimum number of haplotype count.\n\n")

    parser.add_argument("--condition-gene", type=str, help="\nGene name for conditioning.\n\n")
    parser.add_argument("--exclude-composites", action='store_true',
                        help="\nExclude composite amino acids from omnibus test.\n\n")
    parser.add_argument("--output-composites", action='store_true',
                        help="\nOutput composite amino acids in .haplo.txt.\n\n")

    parser.add_argument("--exhaustive", action='store_true',
                        help="\nRun exhaustive testing of AA combinations in a single HLA gene.\n\n")
    parser.add_argument("--exhaustive-aa-pos", type=int,
                        help="\nSpecify a AA pos of the first element(s) of combinations.\n\n")
    parser.add_argument("--exhaustive-min-aa", type=int, default=2,
                        help="\nMinimum number of AA positions to form a combination.\n\n")
    parser.add_argument("--exhaustive-max-aa", type=int, default=2,
                        help="\nMaximum number of AA positions to form a combination.\n\n")
    parser.add_argument("--exhaustive-no-filter", action='store_true',
                        help="\nDon't filter non-increasing AA combinations.\n\n")

    ##### <for Test> #####

    # args = parser.parse_args('--target example/1958BC '
    #                          '--hg 18 '
    #                          '--reference example/HAPMAP_CEU '
    #                          '--hped example/HAPMAP_CEU_HLA.Ggroup.hped '
    #                          '--out tests/MyHLA-TAPAS/myHLA-TAPAS '
    #                          '--pheno example/1958BC.phe'.split(' '))



    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)


    # Main implementation
    HLA_TAPAS(args.target, args.reference, args.hg, args.out, args.hped,
              f_save_intermediates=args.save_intermediates, _dependency=args.dependency,
              b_mem=args.mem, b_nthreads=args.nthreads, b_niterations=args.niterations,
              _tolerated_diff=args.tolerated_diff,
              _pheno=args.pheno, _pheno_name=args.pheno_name, _condition=args.condition, _condition_list=args.condition_list,
              _reference_bim=args.reference_bim, _covar=args.covar, _covar_name=args.covar_name,
              _pop=args.pop, _sex=args.sex, _pcs=args.pcs, _maf_threshold=args.maf_threshold, f_aa_only=args.aa_only,
              f_remove_samples_by_haplo=args.remove_samples_by_haplo,
              f_remove_samples_aa_pattern=args.remove_samples_aa_pattern,
              _min_haplo_count=args.min_haplo_count, _condition_gene=args.condition_gene,
              f_exclude_composites=args.exclude_composites, f_output_composites=args.output_composites,
              f_exhaustive=args.exhaustive, _exhaustive_aa_pos=args.exhaustive_aa_pos,
              _exhaustive_min_aa=args.exhaustive_min_aa, _exhaustive_max_aa=args.exhaustive_max_aa,
              f_exhaustive_no_filter=args.exhaustive_no_filter
              )
