# -*- coding: utf-8 -*-

import os, sys, re
import subprocess
import pandas as pd
from shutil import which
from math import log10



########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "[%s::WARNING]: " % (os.path.basename(__file__))


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


##### Patterns to use #####
p_Omnibus_AA = re.compile(r'AA_([A-Z0-9]+)_(-?\d+)')
p_Omnibus_INS = re.compile(r'INS_AA_([A-Z0-9]+)_(-?\d+)x(-?\d+)')


class HATK_Manhattan():

    def __init__(self, _assoc_result, _out, _hg, **kwargs):


        for item in _assoc_result:
            if not os.path.exists(item):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given association result file can't be founct('{}'). Please check it again." \
                      .format(item))
                sys.exit()

        if not bool(_out):
            print(std_ERROR_MAIN_PROCESS_NAME + "Please check '--out' argument again.\n")
            sys.exit()

        if not bool(_hg):
            print(std_ERROR_MAIN_PROCESS_NAME + "Please check '--hg' argument again.\n")
            sys.exit()


        isOmnibus = list(map(lambda x : x.endswith('.omnibus'), _assoc_result))

        if all(isOmnibus):

            if not bool(kwargs['_HLA']):
                print(std_ERROR_MAIN_PROCESS_NAME + "Which HLA genes to plot must be specified.\n"
                                                    "Please check the '--HLA' argument again.")
                sys.exit()


        self.result = HLAManhattan(_assoc_result, _out, _hg, _topcol=kwargs['_top_color'],
                                _pointsize=kwargs['_point_size'], _yaxis_unit=kwargs['_yaxis_unit'], _HLA=kwargs['_HLA'])


        # self.removeIntermediates(_out)



    def removeIntermediates(self, _out):

        # *.log
        if os.path.exists(_out+'.log'):
            os.system("rm {}".format(_out+'.log'))




def HLAManhattan(_assoc_result, _out, _hg, _topcol="#FF0000", _min_pos="29.60E6", _max_pos="33.2E6",
              _pointsize="15", _yaxis_unit="10", _p_Rscript=which("Rscript"), _HLA=None,
              _p_src='HLAManhattan/src', _p_data='HLAManhattan/data'):



    ########## < Core Variables > ##########

    # Intermediate path.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)

    # Point color
    _topcol = re.escape(_topcol if _topcol.startswith('#') else ("#"+_topcol))


    # Paths
    p_src = _p_src
    p_data = _p_data



    ########## < Dependency and Argument Checking > ##########

    if not isinstance(_assoc_result, list):
        print(std_ERROR_MAIN_PROCESS_NAME + "The argument \"--assoc-result(-ar)\" wansn't given as list. Please check it again.\n")
        sys.exit()


    # Checking whether homogenous association test result files are given.
    f_homogenous = list(map(lambda x : x.endswith(os.path.splitext(_assoc_result[0])[1]), _assoc_result))

    if not all(f_homogenous):
        print(std_ERROR_MAIN_PROCESS_NAME + "There are different types of association test result.\n"
                                            "Please make sure same types of association test results are given.")
        sys.exit()


    # print(std_MAIN_PROCESS_NAME + "Loading \".assoc.logistic\".\n")
    # print("{0} logistic regression results are given.\n".format(len(_assoc_result)))


    l_TOP_LABEL = []
    l_yaxis = []


    if _assoc_result[0].endswith('.assoc.logistic'):

        ##### *.assoc.logistic(PLINK) #####

        for i in range(0, len(_assoc_result)):

            # print("\n[{0}] : {1}\n".format(i, _assoc_result[i]))

            t_ar = pd.read_csv(_assoc_result[i], sep='\s+', header=0, usecols=["SNP", "P"]).dropna().sort_values("P")

            # MARKER_set = t_ar.iloc[:, 0].tolist()
            # print(std_MAIN_PROCESS_NAME + "Marker Labels are \n{0}".format(MARKER_set))

            # print(t_ar.head())

            TOP_LABEL = t_ar.iat[0, 0]
            TOP_LABEL_value = t_ar.iat[0, 1]

            # print("\nTop signal is \"{0}(P-value : {1})\"".format(TOP_LABEL, TOP_LABEL_value))


            temp = -log10(TOP_LABEL_value)

            if temp <= 10:
                _yaxis = 10
            else:

                if temp % 5 == 0:
                    _yaxis = temp
                else:
                    _yaxis = (int(temp/5) + 1)*5


            # print("maximum y-axis is {0}".format(_yaxis))

            l_TOP_LABEL.append(TOP_LABEL)
            l_yaxis.append(str(_yaxis))



        export_ar = ','.join(_assoc_result)
        export_TOP_LABEL = ','.join(l_TOP_LABEL)
        export_yaxis = ','.join(l_yaxis)



        # hg (Human Genome)
        _knownGene = os.path.join(p_data, 'known_genes/known_genes_chr6.hg{0}.txt'.format(_hg))


        # print("\n")


        ########## < Plotting Manhattan > ##########

        # print(std_MAIN_PROCESS_NAME + "Plotting Manhattan.\n")

        command = [_p_Rscript, os.path.join(p_src, "manhattan_HLA_HATK.R"),
                   export_ar, _out,
                   _pointsize, _topcol,
                   _min_pos, _max_pos, export_TOP_LABEL, export_yaxis, _yaxis_unit,
                   _knownGene, p_src]

        command = ' '.join(command)

        try:
            f_log = open(_out+'.log', 'w')
            # print(command)
            subprocess.run(re.split(r'\s+', command), check=True, stdout=f_log, stderr=f_log)

        except subprocess.CalledProcessError:
            # Fail
            print(std_ERROR_MAIN_PROCESS_NAME + "Plotting manhattan failed. Please check its log file('{}')." \
                  .format(_out+'.log'))
            return -1
        else:
            # Succeed
            f_log.close()
            return _out + ".pdf"



    elif _assoc_result[0].endswith('.omnibus'):

        __RESULTS__ = []

        ##### Omnibus Test #####

        # print(std_MAIN_PROCESS_NAME + "Omnibus")

        l_processed_omnibus = []

        for i in range(0, len(_assoc_result)):

            # Header : [Variant, deltaDeviance, deltaDF, N, log10_P, Residues]
            t_ar = pd.read_csv(_assoc_result[i], sep='\s+', header=0) \
                        .dropna()

            # print(t_ar['log10_P'].map(lambda x : 10**x))

            t_ar['P'] = t_ar['log10_P'].map(lambda x : 10**x)


            # Extracting 'HLA' and 'Relative position' of Amino acid marker
            f_AA = t_ar['Variant'].str.match(p_Omnibus_AA)
            df_AA = t_ar['Variant'][f_AA].str.extract(p_Omnibus_AA, expand=True)
            df_AA.columns = ['HLA', 'REL_POS']
            # print("df_AA :\n{}\n".format(df_AA))


            # Extracting 'HLA' and 'Relative position' of INS a.a. marker
            f_INS = t_ar['Variant'].str.match(p_Omnibus_INS)

            if f_INS.any():
                df_INS = t_ar['Variant'][f_INS].str.extract(p_Omnibus_INS, expand=True)
                df_INS = pd.concat([
                    df_INS.iloc[:, 0],
                    df_INS.iloc[:, [1,2]].apply(lambda x : str(x.astype(int).mean()), axis=1)
                ], axis=1)
            else:
                df_INS = pd.DataFrame([]) # Null DataFrame

            df_INS.columns = ['HLA', 'REL_POS']

            # print("df_INS :\n{}\n".format(df_INS))

            df_HLA_RelPos = pd.concat([df_AA, df_INS], axis=0).sort_index()
            df_HLA_RelPos['REL_POS'] = df_HLA_RelPos['REL_POS'].map(lambda x : float(x)) # as float

            t_ar = pd.concat([t_ar, df_HLA_RelPos], axis=1)
            # print(t_ar)


            for i in range(len(_HLA)):

                ### Plotting by each HLA gene.

                t_ar_byHLA = t_ar[t_ar['HLA'] == _HLA[i]]
                # print(t_ar_byHLA.head())

                t_ar_byHLA_2 = t_ar_byHLA[['Variant', 'P', 'REL_POS']].sort_values('P')
                # print(t_ar_byHLA_2.head())


                TOP_LABEL = t_ar_byHLA_2['Variant'].iat[0]
                TOP_LABEL_value = t_ar_byHLA_2['P'].iat[0]

                # print("\nTop signal is \"{0}(P-value : {1})\"".format(TOP_LABEL, TOP_LABEL_value))

                temp = -log10(TOP_LABEL_value)

                if temp <= 10:
                    _yaxis = 10
                else:

                    if temp % 5 == 0:
                        _yaxis = temp
                    else:
                        _yaxis = (int(temp / 5) + 1) * 5

                # print("maximum y-axis is {0}".format(_yaxis))



                # Briefing

                _forManhattn = _out + '.{}.forManhattan.txt'.format(_HLA[i])
                # print(t_ar_byHLA.sort_values('REL_POS'))
                t_ar_byHLA.sort_values('REL_POS').to_csv(_forManhattn, sep='\t', header=True, index=False)




                export_ar = _forManhattn
                export_TOP_LABEL = TOP_LABEL
                export_yaxis = str(_yaxis)

                _min_pos = str(t_ar_byHLA['REL_POS'].min())
                _max_pos = str(t_ar_byHLA['REL_POS'].max())



                ########## < Plotting Manhattan > ##########

                # print(std_MAIN_PROCESS_NAME + "Plotting Manhattan.\n")

                OUT = _out+'.{}'.format(_HLA[i])

                command = [_p_Rscript, os.path.join(p_src, "manhattan_HLA_HATK.Omnibus.R"),
                           export_ar, OUT,
                           _pointsize, _topcol,
                           _min_pos, _max_pos, export_TOP_LABEL, export_yaxis, _yaxis_unit, p_src]

                command = ' '.join(command)

                try:
                    f_log = open(_out + '.log', 'w')
                    # print(command)
                    subprocess.run(re.split(r'\s+', command), check=True, stdout=f_log, stderr=f_log)

                except subprocess.CalledProcessError:
                    # Fail
                    print(std_ERROR_MAIN_PROCESS_NAME + "Plotting Omnibus manhattan failed. Please check its log file('{}')." \
                          .format(_out + '.log'))
                    return -1
                else:
                    # Succeed
                    f_log.close()
                    return _out + ".pdf"


        return __RESULTS__
