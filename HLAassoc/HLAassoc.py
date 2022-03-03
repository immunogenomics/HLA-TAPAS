# -*- coding: utf-8 -*-

import os, sys, re
import subprocess
import pandas as pd
from os.path import exists, join, dirname, basename
from subprocess import DEVNULL, CalledProcessError
from shutil import which

from HLAassoc.src.reverse_map import reverse_map
from HLAassoc.src.unATtrickBGL import unATtrickBGL

std_MAIN_PROCESS_NAME = "[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "[%s::WARNING]: " % (os.path.basename(__file__))



class HLAassoc(object):

    def __init__(self, _MAIN_MENU, _out, _dependency='./dependency', **kwargs):

        """
        ### Function signatures.

        (1) Logistic Regression (`_MAIN_MENU` == 'LOGISTIC')

        __init__(self, _MAIN_MENU, _out, _dependency='./dependency'
            _vcf,
            _reference_bim=None,
            _covar=None, _covar_name=None, _pheno=None, _pheno_name=None,
            _condition=None, _condition_list=None,
            _hped=None, _chped=None, _hat=None,
        )

        (2) Linear Regression (`_MAIN_MENU` == 'LINEAR')

        __init__(self, _MAIN_MENU, _out, _dependency='./dependency'
            _vcf,
            _reference_bim=None,
            _covar=None, _covar_name=None, _pheno=None, _pheno_name=None,
            _condition=None, _condition_list=None,
            _hped=None, _chped=None, _hat=None,
        )

        (2) Omnibus Test (`_MAIN_MENU` == 'OMNIBUS')
        __init__(self, _MAIN_MENU, _out, _dependency='./dependency'
            _vcf=None,
            _file=None,
            _pop=None,
            _phased=None,
            _fam=None,
            _bim=None,
            _pheno=None,
            _covars=None,
            _maf_threshold=0,
            f_aa_only=True,
            _nthreads=1,
            f_remove_samples_by_haplo=False,
            f_remove_samples_aa_pattern=False,
            _condition=None,
            _condition_gene=None,
            f_exclude_composites=False,
            f_output_composites=False,
            f_exhaustive=False,
            _exhaustive_aa_pos=None,
            _exhaustive_min_aa=2,
            _exhaustive_max_aa=2,
            f_exhaustive_no_filter=False
        )

        """

        ### Intermediate path.
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        if bool(dirname(_out)):
            self.out_dirname = dirname(_out)
            os.makedirs(self.out_dirname, exist_ok=True)
        else:
            self.out_dirname = ''




        ### Main Menu
        if _MAIN_MENU == 'LOGISTIC':


            ## Argument Parsing
            _vcf = kwargs['_vcf']
            _reference_bim = kwargs['_reference_bim']
            _covar = kwargs['_covar']
            _covar_name = kwargs['_covar_name']
            _pheno = kwargs['_pheno']
            _pheno_name = kwargs['_pheno_name']
            _condition = kwargs['_condition']
            _condition_list = kwargs['_condition_list']
            _hped = kwargs['_hped']
            _chped = kwargs['_chped']
            # _hat = kwargs['_hat']



            ## Main variables.
            self.assoc_result = None

            self.vcf = None
            self.a1_allele = None
            self.pheno = None
            self.pheno_name = None
            self.covar = None
            self.covar_name = None
            self.condition = None
            self.condition_list = None



            ## Dependency check

            # PLINK
            if exists(join(_dependency, 'plink')):
                self.plink = join(_dependency, 'plink')
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(_dependency))
                sys.exit()



            ## Conversion from VCF to PLINK.
            if not exists(_vcf):
                print(std_ERROR_MAIN_PROCESS_NAME + "Target data can't be found.('{}') Please check '--target/-t' argument again.".format(_vcf))
                sys.exit()
            else:
                if _vcf.endswith('.vcf.gz') or _vcf.endswith('.vcf'):
                    self.vcf = _vcf
                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Given VCF file('{}') doesn't have file extension either '*.vcf.gz' or '*.vcf'. Please check '--target/-t' argument again.".format(_vcf))



            ## Decoding bim file (Whether to use reference bim file or not.)
            # if bool(_reference_bim) and exists(_reference_bim):
            #     self.bim = _reference_bim



            ## Reverse-map HLA marker label.
            # if bool(_hped) and bool(_chped):
            #
            #     if exists(_hped) and exists(_chped):
            #         self.vcf = reverse_map(self.vcf,
            #                                join(self.out_dirname, basename(_out).rstrip('.vcf.gz') + '.rev_mapped.vcf.gz'),
            #                                _hped, _chped)
            #         # print(self.vcf)
            #
            #     else:
            #         if not exists(_hped):
            #             print(std_WARNING_MAIN_PROCESS_NAME + "HPED file can't be found('{}'). Please check '--hped' argument again.\n"
            #                                                   "Skipping Reverse-mapping of HLA marker labels.".format(_hped))
            #         if not exists(_chped):
            #             print(std_WARNING_MAIN_PROCESS_NAME + "CHPED file can't be found('{}'). Please check '--chped' argument again.\n"
            #                                                   "Skipping Reverse-mapping of HLA marker labels.".format(_chped))
            #
            #
            # # elif bool(_hat):
            # #     # if not exists(_hat):
            # #     #     print(std_WARNING_MAIN_PROCESS_NAME + "HAT file can't be found('{}'). Please check '--hped' argument again.\n"
            # #     #                                           "Skipping Reverse-mapping of HLA marker labels.".format(_hat))
            # #
            # #     # Not yet.
            # #     pass



            ## Phenotype file check
            if bool(_pheno):
                if exists(_pheno):
                    # Phenotype file is given.
                    self.pheno = _pheno

                    if bool(_pheno_name):
                        # Which phenotype column name to use is given.
                        self.pheno_name = _pheno_name
                    else:
                        # Phenotype header check.
                        with open(_pheno, 'r') as f_pheno:

                            header = next(f_pheno)
                            l_header = re.split(r'\s+', header.rstrip('\n'))

                            if len(l_header) > 3:
                                # can't determine which phenotype name to use.
                                print(std_ERROR_MAIN_PROCESS_NAME + "can't determined which Phenotype column to use.('{}') "
                                                                    "Please specify which column to use with '--pheno-name' argument." \
                                      .format(_pheno))
                                sys.exit()

                            elif len(l_header) == 3:
                                # Set the only phenotype column name as the one to use.
                                self.pheno_name = l_header[-1]
                                print(std_WARNING_MAIN_PROCESS_NAME + "Using phenotype column '{}' in '{}' file." \
                                      .format(self.pheno_name, self.pheno))

                            else:
                                print(std_ERROR_MAIN_PROCESS_NAME + "Given phenotype file ('{}') has bizarre number of columns. "
                                                                    "Please check '--pheno' argument again." \
                                      .format(_pheno))
                                sys.exit()

                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Phenotype file can't be found('{}'). Please check '--pheno' argument again." \
                          .format(_pheno))

            # else:
            #     # Check whether Phenotype vector is in *.fam file.
            #     # *.fam file check
            #
            #     l_Phe = []
            #
            #     with open(self.fam, 'r') as f_fam:
            #         for line in f_fam:
            #             l = re.split(r'\s+', line.rstrip('\n'))
            #             l_Phe.append(l[5])
            #
            #     isNoPhe = list(map(lambda x : x == '-9', l_Phe))
            #
            #     if all(isNoPhe):
            #         print(std_ERROR_MAIN_PROCESS_NAME + "Phenotype information can't be found in FAM file('{}'). "
            #                                             "Plesae specify Phenotype information with '--pheno' and '--pheno-name' arguments." \
            #               .format(self.fam))
            #         sys.exit()
            # => Checking *.fam file has been deprecated because VCF file is mainly used in HLAassoc.



            ## Covariate file check
            # Not yet



            ## Condition file check
            # Not yet




            ##### Association Test #####
            print(std_MAIN_PROCESS_NAME + "Performing Logistic Regression.")

            self.assoc_result = \
                self.logistic_regression(_out, self.vcf,
                                         _covar=_covar, _covar_name=_covar_name,
                                         _pheno=self.pheno, _pheno_name=self.pheno_name,
                                         _condition=_condition, _condition_list=_condition_list,
                                         f_asBETA=False, _a1_allele=self.a1_allele)

        elif _MAIN_MENU == 'LINEAR':


            ## Argument Parsing
            _vcf = kwargs['_vcf']
            _reference_bim = kwargs['_reference_bim']
            _covar = kwargs['_covar']
            _covar_name = kwargs['_covar_name']
            _pheno = kwargs['_pheno']
            _pheno_name = kwargs['_pheno_name']
            _condition = kwargs['_condition']
            _condition_list = kwargs['_condition_list']
            _hped = kwargs['_hped']
            _chped = kwargs['_chped']
            # _hat = kwargs['_hat']



            ## Main variables.
            self.assoc_result = None

            self.vcf = None
            self.a1_allele = None
            self.pheno = None
            self.pheno_name = None
            self.covar = None
            self.covar_name = None
            self.condition = None
            self.condition_list = None



            ## Dependency check

            # PLINK
            if exists(join(_dependency, 'plink')):
                self.plink = join(_dependency, 'plink')
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(_dependency))
                sys.exit()



            ## Conversion from VCF to PLINK.
            if not exists(_vcf):
                print(std_ERROR_MAIN_PROCESS_NAME + "Target data can't be found.('{}') Please check '--target/-t' argument again.".format(_vcf))
                sys.exit()
            else:
                if _vcf.endswith('.vcf.gz') or _vcf.endswith('.vcf'):
                    self.vcf = _vcf
                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Given VCF file('{}') doesn't have file extension either '*.vcf.gz' or '*.vcf'. Please check '--target/-t' argument again.".format(_vcf))



            ## Decoding bim file (Whether to use reference bim file or not.)
            # if bool(_reference_bim) and exists(_reference_bim):
            #     self.bim = _reference_bim



            ## Reverse-map HLA marker label.
            # if bool(_hped) and bool(_chped):
            #
            #     if exists(_hped) and exists(_chped):
            #         self.vcf = reverse_map(self.vcf,
            #                                join(self.out_dirname, basename(_out).rstrip('.vcf.gz') + '.rev_mapped.vcf.gz'),
            #                                _hped, _chped)
            #         # print(self.vcf)
            #
            #     else:
            #         if not exists(_hped):
            #             print(std_WARNING_MAIN_PROCESS_NAME + "HPED file can't be found('{}'). Please check '--hped' argument again.\n"
            #                                                   "Skipping Reverse-mapping of HLA marker labels.".format(_hped))
            #         if not exists(_chped):
            #             print(std_WARNING_MAIN_PROCESS_NAME + "CHPED file can't be found('{}'). Please check '--chped' argument again.\n"
            #                                                   "Skipping Reverse-mapping of HLA marker labels.".format(_chped))


            ## Phenotype file check
            if bool(_pheno):
                if exists(_pheno):
                    # Phenotype file is given.
                    self.pheno = _pheno

                    if bool(_pheno_name):
                        # Which phenotype column name to use is given.
                        self.pheno_name = _pheno_name
                    else:
                        # Phenotype header check.
                        with open(_pheno, 'r') as f_pheno:

                            header = next(f_pheno)
                            l_header = re.split(r'\s+', header.rstrip('\n'))

                            if len(l_header) > 3:
                                # can't determine which phenotype name to use.
                                print(std_ERROR_MAIN_PROCESS_NAME + "can't determined which Phenotype column to use.('{}') "
                                                                    "Please specify which column to use with '--pheno-name' argument." \
                                      .format(_pheno))
                                sys.exit()

                            elif len(l_header) == 3:
                                # Set the only phenotype column name as the one to use.
                                self.pheno_name = l_header[-1]
                                print(std_WARNING_MAIN_PROCESS_NAME + "Using phenotype column '{}' in '{}' file." \
                                      .format(self.pheno_name, self.pheno))

                            else:
                                print(std_ERROR_MAIN_PROCESS_NAME + "Given phenotype file ('{}') has bizarre number of columns. "
                                                                    "Please check '--pheno' argument again." \
                                      .format(_pheno))
                                sys.exit()

                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Phenotype file can't be found('{}'). Please check '--pheno' argument again." \
                          .format(_pheno))


            ##### Association Test #####
            print(std_MAIN_PROCESS_NAME + "Performing Linear Regression.")

            self.assoc_result = \
                self.linear_regression(_out, self.vcf,
                                         _covar=_covar, _covar_name=_covar_name,
                                         _pheno=self.pheno, _pheno_name=self.pheno_name,
                                         _condition=_condition, _condition_list=_condition_list,
                                          _a1_allele=self.a1_allele)


        elif _MAIN_MENU == 'OMNIBUS':

            ### Argument Parsing

            _vcf = kwargs['_vcf']
            _file = kwargs['_file']
            _pop = kwargs['_pop']
            _phased = kwargs['_phased']
            _fam = kwargs['_fam']
            _bim = kwargs['_bim']
            _pheno = kwargs['_pheno']
            _covars = kwargs['_covars']
            _maf_threshold = kwargs['_maf_threshold']
            f_aa_only = kwargs['f_aa_only']
            _nthreads = kwargs['_nthreads']
            f_remove_samples_by_haplo = kwargs['f_remove_samples_by_haplo']
            f_remove_samples_aa_pattern = kwargs['f_remove_samples_aa_pattern']
            _min_haplo_count=kwargs['_min_haplo_count']
            _condition = kwargs['_condition']
            _condition_gene = kwargs['_condition_gene']
            f_exclude_composites = kwargs['f_exclude_composites']
            f_output_composites = kwargs['f_output_composites']
            f_exhaustive = kwargs['f_exhaustive']
            _exhaustive_aa_pos = kwargs['_exhaustive_aa_pos']
            _exhaustive_min_aa = kwargs['_exhaustive_min_aa']
            _exhaustive_max_aa = kwargs['_exhaustive_max_aa']
            f_exhaustive_no_filter = kwargs['f_exhaustive_no_filter']

            _java_heap_mem = kwargs['_java_heap_mem']



            ### Main Variables

            self.omnibus_result = None

            self.bgl_phased = None
            self.bim = None
            self.fam = None
            self.covars = None
            self.pheno = None
            self.pop = None

            self.Rscript = None


            ### Dependency check.

            # PLINK
            if exists(join(_dependency, 'plink')):
                self.plink = join(_dependency, 'plink')
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(_dependency))
                sys.exit()

            # R
            if exists(which("Rscript")):
                self.Rscript = which("Rscript")
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Please check whether 'Rscript' command is prepared in your system.")
                sys.exit()


            ### Exception Handling

            if bool(_file):
                # Use given common prefix to nominate multiple files at once.
                self.bgl_phased = _file + '.bgl.phased'
                self.bim = _file + '.bim'
                self.fam = _file + '.fam'
                self.covars = _file + '.covs'
                self.pheno = _file + '.pheno'
                self.pop = _file + '.pop'
                #self.sex = _file + '.sex'

            else:

                # *.bim
                if bool(_bim):
                    if exists(_bim):
                        self.bim = _bim
                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Given PLINK bim file('{}') can't be found. "
                                                            "Please check '--bim' argument again." \
                              .format(_bim))
                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "PLINK bim file must be given. "
                                                        "Please specify it with '--bim' argument.")
                    sys.exit()



                # *.fam
                if bool(_fam):
                    if exists(_fam):
                        self.fam = _fam
                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Given PLINK fam file('{}') can't be found. "
                                                            "Please check '--bim' argument again." \
                              .format(_fam))
                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "PLINK fam file must be given. "
                                                        "Please specify it with '--fam' argument.")
                    sys.exit()



                # *.bgl.phased
                if bool(_phased):
                    if exists(_phased):
                        self.bgl_phased = _phased
                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Given Phased BEAGLE file('{}') can't be found. "
                                                            "Please check '--phased' argument again." \
                              .format(_phased))
                        sys.exit()

                elif not bool(_phased) and bool(_vcf):
                    if exists(_vcf):
                        # Conversion from VCF to BEAGLE file. (vcf2beagle.jar)
                        print(std_MAIN_PROCESS_NAME + "Phased BEAGLE file will be generated from given VCF file('{}')." \
                              .format(_vcf))

                        if not exists(join(_dependency, 'vcf2beagle.jar')):
                            print(std_ERROR_MAIN_PROCESS_NAME +
                                  "'vcf2beagle.jar' has be prepared in 'dependency/' folder to convert given VCF file"
                                  "('{}') to Phased BEAGLE file.".format(_vcf))
                            sys.exit()

                        self.bgl_phased = self.VCF2BEAGLE(_vcf, join(self.out_dirname, re.sub(r'\.vcf\.gz$', '', basename(_vcf))),
                                                          self.bim, self.fam,
                                                          join(_dependency, 'vcf2beagle.jar'), _java_heap_mem)

                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Given VCF file('{}') can't be found. "
                                                            "Please check '--vcf' argument again." \
                              .format(_vcf))
                        sys.exit()

                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "Target file for Omnibus Test must be given. "
                                                        "Please specify it with '--file', '--vcf' or '--phased' argument.")
                    sys.exit()



                # *.covs covariate file
                if bool(_covars):
                    if exists(_covars):
                        self.covars = _covars
                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Given covariate information file('{}') can't be found. "
                                                            "Please check '--covars' argument again." \
                              .format(_covars))


                elif not bool(_covars) and bool(_vcf):

                    print(std_MAIN_PROCESS_NAME + "Top 10 PCs will be generated from given VCF file('{}').".format(_vcf))

                    if exists(_vcf):

                        self.covars = self.getTOP10PCs(_vcf, _out, self.plink)

                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Given VCF file('{}') to generate top 10 PCs can't be found." \
                              .format(_vcf))
                        sys.exit()


                else:
                    print(std_ERROR_MAIN_PROCESS_NAME + "covariate information file must be given. "
                                                        "Please specify it with '--covars' argument.")
                    sys.exit()



                # *.pheno
                if bool(_pheno):
                    if exists(_pheno):
                        self.pheno = _pheno
                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Given Phenotype file('{}') can't be found. "
                                                          "Please check '--pheno' argument again." \
                            .format(_pheno))

                else:

                    if self.hasPHENOTYPEinFAM(self.fam):

                        print(std_WARNING_MAIN_PROCESS_NAME +
                              "Phenotype information in given fam file('{}') will be used.".format(self.fam))

                        df_temp_phe = pd.read_csv(self.fam, sep='\s+', header=None, dtype=str,
                                                  names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe'],
                                                  usecols=['FID', 'IID', 'Phe'])


                        if df_temp_phe['Phe'].map(lambda x: (x == '1' or x == '2')).all():

                            df_temp_phe \
                                .pipe(self.ZeroOnePhenotype) \
                                .to_csv(_out+'.pheno', sep=' ', header=False, index=False)

                        else:
                            df_temp_phe.to_csv(_out+'.pheno', sep=' ', header=False, index=False)


                        self.pheno = _out+'.pheno'


                    else:
                        print(std_ERROR_MAIN_PROCESS_NAME + "Phenotype file must be given. "
                                                            "Please specify it with '--pheno' argument.")
                        sys.exit()


                # *.pop
                # if bool(_pop):
                #    if exists(_pop):
                #        self.pop = _pop
                #    else:
                #        print(std_ERROR_MAIN_PROCESS_NAME + "Given Population information file('{}') can't be found. "
                #                                            "Please check '--pop' argument again." \
                #            .format(_pop))
		#
                # else:
                #    print(std_WARNING_MAIN_PROCESS_NAME +
                #         "All samples will be assumed to be originated from same population.")

                #    pd.read_csv(self.fam, sep='\s+', header=None, dtype=str,
                #                names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe']) \
                #        .pipe(self.getDummyPOP) \
                #        .to_csv(_out+'.pop', sep=' ', header=False, index=False)

                 #   self.pop = _out+'.pop'

                    # print(std_ERROR_MAIN_PROCESS_NAME + "Population information file must be given. "
                    #                                     "Please specify it with '--pop' argument.")
                    # sys.exit()




            ##### Omnibus Test #####

            print(std_MAIN_PROCESS_NAME + "Performing Omnibus Test.")
            self.omnibus_result = self.Omnibus_Test(
                _out, self.pop, self.bgl_phased, self.fam, self.bim, self.pheno,  self.covars,
                _maf_threshold, f_aa_only, _nthreads, f_remove_samples_by_haplo, f_remove_samples_aa_pattern,
                _min_haplo_count,
                _condition, _condition_gene, f_exclude_composites, f_output_composites, f_exhaustive,
                _exhaustive_aa_pos, _exhaustive_min_aa, _exhaustive_max_aa, f_exhaustive_no_filter
            )



        else:

            pass

        # End __init__()





    ##### Class methods  #####

    def logistic_regression(self, _out, _vcf,
                            _covar=None, _covar_name=None, _pheno=None, _pheno_name=None,
                            _condition=None, _condition_list=None, _a1_allele=None,
                            f_asBETA=False, _ci=0.95):


        command = [self.plink, '--vcf {}'.format(_vcf), '--double-id',
                   '--allow-no-sex', '--ci {}'.format(_ci), '--out {}'.format(_out)]


        # as Beta vs. OR ?
        if f_asBETA:
            command.extend(['--logistic', 'hide-covar', 'beta'])
        else:
            command.extend(['--logistic', 'hide-covar'])

        # a1_allele check
        if bool(_a1_allele):
            command.append('--a1-allele {}'.format(_a1_allele))

        # Phenotype check
        if bool(_pheno):
            command.append('--pheno {}'.format(_pheno))

            if bool(_pheno_name):
                command.append('--pheno-name {}'.format(_pheno_name))

        # Covariate check
        if bool(_covar):
            command.append('--covar {}'.format(_covar))

            if bool(_covar_name):
                command.append('--covar-name {}'.format(_covar_name))

        # Condition check
        if bool(_condition):
            command.append('--condition {}'.format(_condition))
        elif bool(_condition_list):
            command.append('--condition-list {}'.format(_condition_list))


        try:
            # print(command)
            subprocess.run(re.split(r'\s+', ' '.join(command)), check=True, stdout=DEVNULL, stderr=DEVNULL)
        except CalledProcessError:
            # Fail
            print(std_ERROR_MAIN_PROCESS_NAME + "Association Test Failed. Check plink log file for more details")
            sys.exit()
        else:
            # Succeed

            if exists(_out+'.nosex'):
                os.system('rm {}'.format(_out+'.nosex'))

            return _out+'.assoc.logistic'



    def linear_regression(self, _out, _vcf,
                            _covar=None, _covar_name=None, _pheno=None, _pheno_name=None,
                            _condition=None, _condition_list=None, _a1_allele=None, _ci=0.95):


        command = [self.plink, '--vcf {}'.format(_vcf),"--double-id",
                   '--allow-no-sex', '--ci {}'.format(_ci), '--out {}'.format(_out)]


        # perform linear regression
        command.extend(['--linear', 'hide-covar', 'beta'])

        # a1_allele check
        if bool(_a1_allele):
            command.append('--a1-allele {}'.format(_a1_allele))

        # Phenotype check
        if bool(_pheno):
            command.append('--pheno {}'.format(_pheno))

            if bool(_pheno_name):
                command.append('--pheno-name {}'.format(_pheno_name))

        # Covariate check
        if bool(_covar):
            command.append('--covar {}'.format(_covar))

            if bool(_covar_name):
                command.append('--covar-name {}'.format(_covar_name))

        # Condition check
        if bool(_condition):
            command.append('--condition {}'.format(_condition))
        elif bool(_condition_list):
            command.append('--condition-list {}'.format(_condition_list))


        try:
            # print(command)
            subprocess.run(re.split(r'\s+', ' '.join(command)), check=True, stdout=DEVNULL, stderr=DEVNULL)
        except CalledProcessError:
            # Fail
            print(std_ERROR_MAIN_PROCESS_NAME + "Association Test Failed.")
            sys.exit()
        else:
            # Succeed

            if exists(_out+'.nosex'):
                os.system('rm {}'.format(_out+'.nosex'))

            return _out+'.assoc.linear'

    def Omnibus_Test(self, _out, _pop, _bgl_phased, _fam, _bim, _pheno, _covars,
                     _maf_threshold=0,
                     f_aa_only=True,
                     _nthreads=1,
                     f_remove_samples_by_haplo=False,
                     f_remove_samples_aa_pattern=False,
                     _min_haplo_count=10,
                     _condition=None,
                     _condition_gene=None,
                     f_exclude_composites=True,
                     f_output_composites=False,
                     f_exhaustive=False,
                     _exhaustive_aa_pos=None,
                     _exhaustive_min_aa=2,
                     _exhaustive_max_aa=2,
                     f_exhaustive_no_filter=False,
                     _p_script='HLAassoc/src/run_omnibus_test_LOGISTIC.R'):


        necessary = "--out {} --pop {} --phased {} --fam {} --bim {} --pheno {}  --covars {}" \
                    .format(_out, _pop, _bgl_phased, _fam, _bim, _pheno, _covars)

        default = "--maf-threshold {} --n-threads {} --exhaustive-min-aa {} --exhaustive-max-aa {} --min-haplo-count {}" \
                    .format(_maf_threshold, _nthreads, _exhaustive_min_aa, _exhaustive_max_aa, _min_haplo_count)


        optional = "--exhaustive" if f_exhaustive else "--omnibus"

        if f_aa_only:
            optional = optional + " --aa-only"
        if f_remove_samples_by_haplo:
            optional = optional + " --remove-samples-by-haplo"
        if f_remove_samples_aa_pattern:
            optional = optional + " --remove-samples-aa-pattern"
        if bool(_condition):
            optional = optional + " --condition {}".format(_condition)
        if bool(_condition_gene):
            optional = optional + " --condition-gene {}".format(_condition_gene)
        if f_exclude_composites:
            optional = optional + " --exclude-composites"
        if f_output_composites:
            optional = optional + " --output-composites"
        if bool(_exhaustive_aa_pos):
            optional = optional + " --exhaustive-aa-pos {}".format(_exhaustive_aa_pos)
        if f_exhaustive_no_filter:
            optional = optional + " --exhaustive-no-filter"


        command = ' '.join([self.Rscript, _p_script, necessary, default, optional])

        try:
            f_log = open(_out+'.OMlog', 'w')
            # print(command)
            subprocess.run(re.split(r'\s+', command), check=True, stdout=f_log, stderr=f_log)

        except subprocess.CalledProcessError:
            # Fail
            print(std_ERROR_MAIN_PROCESS_NAME + "Omnibus Test failed. See the log file('{}')." \
                  .format(_out+'.OMlog'))
        else:
            # Succeed
            f_log.close()
            return _out




    def Make_a1_allele(self, _bim, _a1_allele):

        command = 'awk \'{print $2"\t"$5}\' %s > %s' % (_bim, _a1_allele)
        # print(command)
        r = os.system(command)

        if r == 0:
            return _a1_allele
        else:
            print(std_WARNING_MAIN_PROCESS_NAME + "Generating *.a1_allele file failed.")




    def VCF2PLINK(self, _target_vcf, _out):

        command = ' '.join([self.plink,
                            '--vcf {}'.format(_target_vcf),
                            '--make-bed',
                            '--out {}'.format(_out),
                            '--allow-no-sex'])

        try:
            # print(command)
            subprocess.run(re.split(r'\s+', command), check=True, stdout=DEVNULL, stderr=DEVNULL)
        except CalledProcessError:
            # Fail
            print(std_ERROR_MAIN_PROCESS_NAME + "Conversion from VCF to PLINK failed.")
            sys.exit()
        else:
            # Succeed
            return _out




    def VCF2BEAGLE(self, _vcf, _out, _ref_bim, _tar_fam,
                   _vcf2beagle, _mem):

        if _vcf.endswith('.vcf.gz'):
            command1 = 'gunzip -c {}'.format(_vcf)
        else:
            command1 = 'cat {}'.format(_vcf)


        command2 = 'java -Xmx{} -jar {} 0 {}'.format(_mem, _vcf2beagle, _out)

        command = ' | '.join([command1, command2])
        # print(command)
        r = os.system(command)

        if r != 0:
            print(std_ERROR_MAIN_PROCESS_NAME + "Converting given VCF file('{}') to Phased BEAGLE file failed." \
                  .format(_vcf))
            sys.exit()


        raw_bgl = _out+'.bgl.gz'

        # Applying un-ATtrick.
        BGL = unATtrickBGL(raw_bgl, _ref_bim, _out, _tar_fam)


        # Removal
        os.system('rm {}'.format(raw_bgl))
        os.system('rm {}'.format(_out+'.int'))
        os.system('rm {}'.format(_out+'.markers'))

        return BGL



    # remove this function since pcs should not be calculated based on MHC region only. Yang (March 03. 2022)
    # def getTOP10PCs(self, _vcf, _out, _plink):
    #
    #    command = '{} --pca 10 header tabs --allow-no-sex --vcf {} --out {}'.format(_plink, _vcf, _out)

     #   try:
     #       # print(command)
      #      subprocess.run(re.split(r'\s+', command), check=True, stdout=DEVNULL, stderr=DEVNULL)
      #  except CalledProcessError:
     #       # Fail
     #       print(std_ERROR_MAIN_PROCESS_NAME + "Generating Top 10 PCs failed. See PLINK log file('{}')." \
     #             .format(_out+'.log'))
     #       sys.exit()
     #   else:
      #      # Succeed
       #     os.system('mv {} {}'.format(_out+'.eigenvec', _out+'.covs'))
        #    os.system('rm {}'.format(_out+'.eigenval'))
        #    os.system('rm {}'.format(_out+'.log'))
         #   if exists(_out+'.nosex'):
          #      os.system('rm {}'.format(_out+'.nosex'))

           # return _out+'.covs'




    def hasPHENOTYPEinFAM(self, _fam):

        df_fam = pd.read_csv(_fam, sep='\s+', header=None, dtype=str, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe'])

        f_NA1 = df_fam['Phe'] == '-9'
        f_NA2 = df_fam['Phe'].isna()

        f_RETURN = f_NA1 | f_NA2

        return (not f_RETURN.all())




    def hasSEXinFAM(self, _fam):

        df_fam = pd.read_csv(_fam, sep='\s+', header=None, dtype=str, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe'])

        f_NA1 = df_fam['Sex'] == '-9'
        f_NA2 = df_fam['Sex'] == '0'
        f_NA3 = df_fam['Sex'].isna()

        f_RETURN = f_NA1 | f_NA2 | f_NA3

        return (not f_RETURN.all())




    def hasFIDHeader(self, _file):

        with open(_file, 'r') as f:

            line_1st = f.readline()
            return line_1st.startswith('FID')




    def ZeroOnePhenotype(self, _df_phe):

        return pd.concat([_df_phe.iloc[:, [0,1]], (_df_phe.iloc[:, 2].astype(int) - 1)], axis=1)




    def getDummyPOP(self, _df_fam):

        return pd.concat([_df_fam.iloc[:, [0, 1]],
                          pd.Series(['dummy_pop']*_df_fam.shape[0], name='POP', index=_df_fam.index)], axis=1)
