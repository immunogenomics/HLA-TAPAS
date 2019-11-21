# -*- coding: utf-8 -*-

# (2017/11/21) Initiated recoding by Wanson Choi

# (2017.12.1) I replaced using '+=' operator to using '.append()' function. I found that '+=' operator is too slow because it entails copying whole contents of previous one.


import os, sys
import pandas as pd
import argparse, textwrap


def HLAtoSequences(_p_ped, _dictionary, _type, _out):

    ########## < Core Variables > ##########

    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

    # isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}

    HLA_DICTIONARY = pd.DataFrame()
    HLA_DICTIONARY_dict = {}
    ALLELES_SEQ_LENGTH = {}

    # Module Name for stdout
    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME + "Init.")


    ########## < Argument checking > ##########

    # (1) ped file existence
    if not os.path.isfile(_p_ped):
        print(std_MAIN_PROCESS_NAME + "Given ped file doen't exist. Please check it againg.\n")
        sys.exit()

    # (2) HLA DICTIONARY file
    if not os.path.isfile(_dictionary):
        print(std_MAIN_PROCESS_NAME + "Given dictionary file doen't exist. Please check it againg.\n")
        sys.exit()

    # (3) Chekcing `_type`
    if not (_type == "AA" or _type == "SNPS"):
        print(std_MAIN_PROCESS_NAME + "Given value for argument `_type` has wrong value. Please check it againg.\n")
        sys.exit()


    ########## < Control Flags > ##########

    LOADING_DICTIONARY = 1
    LOADING_COATED_PED = 1
    BRINGING_SEQUENCES = 1
    EXPORTING_OUTPUT_PED = 1



    if LOADING_DICTIONARY:

        ########## <1. Dictionary Preparation> ##########

        print("\n[1] Loading Dictionary Data.\n")

        if os.path.isfile(_dictionary):
            HLA_DICTIONARY = pd.read_table(_dictionary, sep='\t', header=None, names=["Alleles", "Seqs"], index_col=0)
        else:
            print(std_MAIN_PROCESS_NAME + "Given Dictionary file doesn't exit!\n")
            sys.exit()

        # (2018. 7. 13.) deprecated - going back to use HLA gene captioned way.
        # # Processing HLA gene caption parts(splited by '*')
        # df_temp = pd.DataFrame(HLA_DICTIONARY.loc[:, "Alleles"].apply(lambda x : x.split('*')).tolist(), columns=["HLA", "Alleles"])
        #
        # HLA_DICTIONARY = pd.concat([df_temp, HLA_DICTIONARY.loc[:, "Seqs"]], axis=1).set_index('HLA')
        print(HLA_DICTIONARY.head())


        ### Dividing `HLA_DICTIONARY` in advance.

        # For performance efficiency, `HLA_DICTIONARY` will be divded by HLA gene type in advance.
        HLA_DICTIONARY_dict = {HLA_names[i]: HLA_DICTIONARY.filter(regex= ''.join(["^", HLA_names[i], "\*"]), axis=0) for i in range(0, len(HLA_names))}
        # HLA_DICTIONARY_dict = {HLA_names[i]: HLA_DICTIONARY.loc[HLA_names[i], :].reset_index(drop=True).set_index('Alleles') for i in range(0, len(HLA_names))}

        for i in range(0, len(HLA_names)):
            print("\nSequence Information of %s" % HLA_names[i])
            print(HLA_DICTIONARY_dict[HLA_names[i]].head())

        ALLELES_SEQ_LENGTH = {HLA_names[i] : len(HLA_DICTIONARY_dict[HLA_names[i]].iat[0, 0]) for i in range(0, len(HLA_names))}

        for i in range(0, len(HLA_names)):
            print("\nSequence Length of %s" % HLA_names[i])
            print(ALLELES_SEQ_LENGTH[HLA_names[i]])



    if LOADING_COATED_PED:

        ########## <2. Loading Coated PED(Input PED) file> ##########

        print("\n[2] Loading Input PED file.")

        INPUT_PED = pd.read_table(_p_ped, sep='\t', header=None, index_col=[0, 1, 2, 3, 4, 5], dtype=str)
        INPUT_PED.columns = pd.Index([name + '_' + str(j + 1) for name in HLA_names for j in range(0, 2)])

        print(INPUT_PED.head())



    if BRINGING_SEQUENCES:

        ########## <3. Bringing Sequences> ##########

        print("\n[3]: Transforming Allele names to Sequences.")

        l_FOR_OUTPUT = []
        l_FOR_OUTPUT_test = []

        for i in range(0, len(HLA_names)):

            curr_dict = HLA_DICTIONARY_dict[HLA_names[i]]

            df_temp = INPUT_PED.filter(regex='_'.join([HLA_names[i], '\d{1}']), axis=1)
            # print(df_temp.head())

            df_temp = df_temp.applymap(lambda x : BringSequence(x, curr_dict) if x != "0" else x)
            # print(df_temp.head())

            l_FOR_OUTPUT_test.append(df_temp)

            # print("\n===============\n")

            # Now, we need to iterate on the number of rows of this DataFrame

            COLUMNS_BY_HLA = []
            # COLUMNS_BY_HLA_test = []

            for j in range(0, len(df_temp)):

                if df_temp.iat[j, 0] != '0' and df_temp.iat[j, 1] != '0':

                    # (Case1) Most normal case - wehn allele_name is found as pair.
                    # ex) allele1 : A*25:01:01  /  allele2 : A*03:01:01:01

                    # seq1 = df_temp.iat[j, 0] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 0][::-1]
                    # seq2 = df_temp.iat[j, 1] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 1][::-1]

                    # (2018. 3. 9) 다시 여기서 reverse안시키는 걸로 바꿈
                    seq1 = df_temp.iat[j, 0]
                    seq2 = df_temp.iat[j, 1]

                    PAIRED = [value for item in zip(seq1, seq2) for value in item]

                else:

                    # (Case2) when not found as a pair of alleles, but as a only single allele, => 0 is given
                    # (0, 0) will be compensated as length of `HLA_seq`, Due to here, I need to prepared `len(HLA_seq)` beforehand.

                    seq1 = ['0' for z in range(0, ALLELES_SEQ_LENGTH[HLA_names[i]])]

                    PAIRED = [value for item in zip(seq1, seq1) for value in item]

                COLUMNS_BY_HLA.append(PAIRED)
                # COLUMNS_BY_HLA_test.append(''.join(PAIRED))


            l_FOR_OUTPUT.append(pd.DataFrame(COLUMNS_BY_HLA))
            # l_FOR_OUTPUT_test.append(pd.Series(COLUMNS_BY_HLA_test))

            # End of interation. Don't get confused.



    if EXPORTING_OUTPUT_PED:

        ########## <4. Exporting OUTPUT PED file> ##########

        print("\n[4]: Exporting OUTPUT PED file.")

        df_OUTPUT = pd.concat(l_FOR_OUTPUT, axis=1)
        df_OUTPUT.index = INPUT_PED.index
        df_OUTPUT.columns = pd.Index(range(0, df_OUTPUT.shape[1]))

        print(df_OUTPUT.head())



        ### Final Output ped file.
        if _type == 'AA':
            df_OUTPUT.to_csv(_out + '.AA.ped', sep='\t', header=False, index=True)
        elif _type == 'SNPS':
            df_OUTPUT.to_csv(_out + '.SNPS.ped', sep='\t', header=False, index=True)


    return 0

def BringSequence(_single_allele, _dict):

    try:
        Seq = _dict.loc[_single_allele, "Seqs"]
    except KeyError:
        Seq = "0"

    return Seq

if __name__ == '__main__':

    # Introduced argparse 도입(2017.11.24)

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

        Original Author: Sherman Jia, 2012

        HLAtoSequences.py
        - This script Converts HLA alleles (in .ped file format) to amino acid or DNA sequences

        Input file should contain: FID, IID, pID, mID, sex, pheno, HLA-A (2), B (2), C (2), 
        DPA1 (2), DPB1 (2), DQA1 (2), DQB1 (2), DRB1 (2) ... Broad Order

    #########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"
    # parser._optionals.description = "- Necessary main options.\n"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-dict", help="\nHLA dictonary file name(ex. 'HLA_DICTIONARY_AA.txt')\n\n", required=True)
    parser.add_argument("-type", help="\nAA(for Amino Acid) or SNP(for SNPs)\n\n", choices=["AA", "SNPS"], required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)




    ##### <for Test> #####

    # (2018.2.9)
    # args = parser.parse_args(["./HAPMAP_CEU_HLA_switched.ped", "./HLA_DICTIONARY_AA_hg19.txt", "AA", "--OUTPUT", "BROAD_ORDER"])


    # (2018.2.26)
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./HLA_DICTIONARY_AA.hg19.imgt370.txt", "AA", "--OUTPUT", "TEST_0329_HLAtoSeq"])

    # (2018.3.8)
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./HLA_DICTIONARY_SNPS.hg19.imgt370.txt", "SNPS", "--OUTPUT", "BROAD_ORDER"])

    # (2018. 7. 12.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370"])


    ##### <for Publication> #####

    args = parser.parse_args()


    print(args)


    # Implementing Main Function
    HLAtoSequences(args.ped, args.dict, args.type, args.o)


