# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names2 = [0, 2, 1, 7, 5, 6, 3, 4]  # Read in the order of `HLA_names`, Write in the order of `HLA_names2` (to be compatible with old version of Dictionary).

isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}


def HLAtoSequences(_chped, _dictionary, _type, _out, __previous_version=False, __asLump=False):



    ### Intermediate path.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)):
        INTERMEDIATE_PATH = os.path.dirname(_out)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        INTERMEDIATE_PATH = "./"


    ########## < Argument checking > ##########

    # (1) ped file existence
    if not os.path.isfile(_chped):
        print(std_MAIN_PROCESS_NAME + "Given ped file doen't exist. Please check it again.\n")
        sys.exit()

    # (2) HLA DICTIONARY file
    if not os.path.isfile(_dictionary):
        print(std_MAIN_PROCESS_NAME + "Given dictionary file doen't exist. Please check it againg.\n")
        sys.exit()

    # (3) Chekcing `_type`
    if not (_type == "AA" or _type == "SNPS"):
        print(std_MAIN_PROCESS_NAME + "Given value for argument `_type` has wrong value. Please check it again.\n")
        sys.exit()



    ##### < Core Variables > #####

    HLA_DICTIONARY_byHLA = {HLA_names[i]: {} for i in range(0, len(HLA_names))}  # Initialization.
    HLA_SEQ_LENGTH = {HLA_names[i]: -1 for i in range(0, len(HLA_names))}

    HLA_INS_byHLA = {HLA_names[i]: {} for i in range(0, len(HLA_names))}  # Initialization.
    haveInsertion= {HLA_names[i]: -1 for i in range(0, len(HLA_names))}


    if __previous_version:

        ##### < [1] Loading HLA Dictionary > #####

        with open(_dictionary, "r") as f_dictionary:

            count = 0

            for line in f_dictionary:

                t_line = re.split(r'\s+', line.rstrip('\n'))
                # ex1) (AA) ['B*58:01:01', 'MLVMAPRTVLLLLSAALALTETWAG...', 'x'] (len == 3)
                # ex2) (SNPS) ['B*58:01:01', 'CTAGTCCTGCTTCAGGGTCCGGGGCCCG...'] (len == 2)
                # print(t_line)

                for i in range(0, len(HLA_names)):

                    if re.match(r'{}:'.format(HLA_names[i]), t_line[0]):

                        HLA_DICTIONARY_byHLA[HLA_names[i]][t_line[0]] = t_line[1] # Sequence information.

                        if HLA_SEQ_LENGTH[HLA_names[i]] == -1:
                            HLA_SEQ_LENGTH[HLA_names[i]] = len(t_line[1])


                        if _type == "AA":

                            HLA_INS_byHLA[HLA_names[i]][t_line[0]] = t_line[2] # Insertion information.

                            if haveInsertion[HLA_names[i]] == -1:
                                haveInsertion[HLA_names[i]] = bool(t_line[2])

                        break  # If a line of given dictionary belongs to either HLA, then checking whether it belongs to other HLAs is useless.

                count += 1
                # if count > 5 : break

        # # Result check
        # for i in range(0, len(HLA_names)):
        #     print("\n{} :\n".format(HLA_names[i]))
        #     for k, v in HLA_DICTIONARY_byHLA[HLA_names[i]].items():
        #         print("{} : {}".format(k, v))
        #
        # for k, v in HLA_SEQ_LENGTH.items():
        #     print("The length of HLA-{} : {}".format(k, v))
        #
        # print("Insertion check : {}".format(haveInsertion))


        ##### < [2] Transforming each HLA alleles to corresponding sequences > #####

        with open(_out + ".{}.ped".format(_type), 'w') as f_output:
            f_output.writelines(GenerateLines(_chped, _type, HLA_DICTIONARY_byHLA, HLA_SEQ_LENGTH, HLA_INS_byHLA, haveInsertion,
                                              __asLump=__asLump))


    else:

        ##### < [1] Loading HLA Dictionary > #####

        with open(_dictionary, "r") as f_dictionary:

            count = 0

            for line in f_dictionary:

                t_line = re.split(r'\s+', line.rstrip('\n'))
                # ex1) (AA) ['B*58:01:01', 'MLVMAPRTVLLLLSAALALTETWAG...'] (len == 2)
                # ex2) (SNPS) ['B*58:01:01', 'CTAGTCCTGCTTCAGGGTCCGGGGCCCG...'] (len == 2)
                # print(t_line)

                for i in range(0, len(HLA_names)):

                    if re.match(r'{}\*'.format(HLA_names[i]), t_line[0]):

                        HLA_DICTIONARY_byHLA[HLA_names[i]][t_line[0]] = t_line[1] # Sequence information.

                        if HLA_SEQ_LENGTH[HLA_names[i]] == -1:
                            HLA_SEQ_LENGTH[HLA_names[i]] = len(t_line[1])

                        break  # If a line of given dictionary belongs to either HLA, then checking whether it belongs to other HLAs is useless.

                count += 1
                # if count > 5 : break


        # # Result check
        # for i in range(0, len(HLA_names)):
        #     print("\n{} :".format(HLA_names[i]))
        #
        #     idx = 0
        #     for k, v in HLA_DICTIONARY_byHLA[HLA_names[i]].items():
        #         print("{} : {}".format(k, v))
        #
        #         idx += 1
        #         if idx > 10 : break
        #
        # for k, v in HLA_SEQ_LENGTH.items():
        #     print("The length of HLA-{} : {}".format(k, v))

        # The insertion will be precessed in the dictionary generation in advance.
        # print("Insertion check : {}".format(haveInsertion))


        ##### < [2] Transforming each HLA alleles to corresponding sequences > #####

        with open(_out + ".{}.ped".format(_type), 'w') as f_output:
            f_output.writelines(GenerateLines(_chped, _type, HLA_DICTIONARY_byHLA, HLA_SEQ_LENGTH, HLA_INS_byHLA, haveInsertion,
                                              __previous_version=__previous_version, __asLump=__asLump))




def GenerateLines(_chped, _type, _dict_seq, _seq_length, _dict_ins, _haveIns,
                  __previous_version=True, __asLump=False):

    with open(_chped, "r") as f:

        for line in f:
            t_line = re.split(r'\s+', line.rstrip('\n'))
            # print(t_line)

            """
            [0,1,2,3,4,5] := ped file information
            [6,7] := HLA-A,
            [8,9] := HLA-B,
            ...,
            [20, 21] := HLA-DRB1
            """

            t_iterator = HLA_names2 if __previous_version else range(0, len(HLA_names))

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_part__ = '\t'.join([
                BringSequence2(t_line[(2 * i + 6)], t_line[(2 * i + 7)], _type, HLA_names[i],
                               _dict_seq[HLA_names[i]], _seq_length[HLA_names[i]],
                               _dict_ins[HLA_names[i]], _haveIns,
                               __previous_version=__previous_version, __asLump=__asLump) for i in t_iterator
            ])

            # mem_p2 = process.memory_info().rss / 1024 ** 2
            # print("{}(Mb)".format(mem_p2 - mem_p1))

            yield '\t'.join([__ped_info__, __genomic_part__]) + "\n"



def BringSequence2(_HLA_allele1, _HLA_allele2, _type, _hla,
                   _dict_seq, _seq_length,
                   _dict_ins, _haveIns,
                   __previous_version = True,
                   __asLump=False):


    if __previous_version:

        # print("al1 : {}\nal2 : {}".format(_HLA_allele1, _HLA_allele2))

        # Finding the corresponding sequence and insertion marker of `_HLA_allele1`.
        try:
            Seq1 = _dict_seq[_HLA_allele1]

            if _type == "AA" and isREVERSE[_hla]:
                Seq1 = Seq1[::-1]

        except KeyError:
            Seq1 = "-1" # fail


        # Same job for `_HLA_allele2`.
        try:
            Seq2 = _dict_seq[_HLA_allele2]

            if _type == "AA" and isREVERSE[_hla]:
                Seq2 = Seq2[::-1]

        except KeyError:
            Seq2 = "-1" # fail


        # print("Corresponding Seqs : \n{}\n{}".format(Seq1, Seq2))

        ### Main sequence information processing
        if not __asLump:

            if Seq1 != "-1" and Seq2 != "-1":

                # Only when both HLA alleles can get the corresponding HLA sequence information.

                l_temp = []

                for i in range(0, len(Seq1)):
                    l_temp.append(Seq1[i])
                    l_temp.append(Seq2[i])

                # Reversing
                __return__ = '\t'.join(l_temp)

            else:
                __return__ = '\t'.join(["0" for z in range(0, 2 * _seq_length)])



            ### Insertion part processing.
            if _type == "AA" and _haveIns[_hla]:

                # One important assumption : "Every insertion part is just one, single character.(ex 'x', 'P', 'A')"
                try:
                    ins1 = _dict_ins[_HLA_allele1]
                except KeyError:
                    ins1 = ""

                try:
                    ins2 = _dict_ins[_HLA_allele2]
                except KeyError:
                    ins2 = ""

                if bool(ins1) and bool(ins2):
                    __insertions__ = '\t'.join([ins1, ins2])
                else:
                    __insertions__ = '\t'.join(["0", "0"])


                __return__ = '\t'.join([__return__, __insertions__])

            return __return__

        else:

            # As each Strings.

            t_Seq1 = ""
            t_Seq2 = ""

            if Seq1 != "-1" and Seq2 != "-1":
                t_Seq1 = Seq1
                t_Seq2 = Seq2
            else:
                t_Seq1 = ["0" for z in range(0, _seq_length)]
                t_Seq2 = ["0" for z in range(0, _seq_length)]

            ### Insertion part processing.
            if _type == "AA" and _haveIns[_hla]:

                # One important assumption : "Every insertion part is just one, single character.(ex 'x', 'P', 'A')"
                try:
                    ins1 = _dict_ins[_HLA_allele1]
                except KeyError:
                    ins1 = ""

                try:
                    ins2 = _dict_ins[_HLA_allele2]
                except KeyError:
                    ins2 = ""


                if bool(ins1):
                    t_Seq1= ''.join([t_Seq1, ins1])

                if bool(ins2):
                    t_Seq2 = ''.join([t_Seq2, ins2])

            return str([[_HLA_allele1, t_Seq1], [_HLA_allele2, t_Seq2]])



    else:

        ### Dealing with generalized 4-field HLA alleles.

        # `_HLA_allele1`.
        try:
            Seq1 = _dict_seq[_HLA_allele1]
        except KeyError:
            Seq1 = "-1"  # fail

        # `_HLA_allele2`.
        try:
            Seq2 = _dict_seq[_HLA_allele2]
        except KeyError:
            Seq2 = "-1"  # fail

        # No reversing HLA sequences.


        if not __asLump:

            ### Main sequence information processing
            if Seq1 != "-1" and Seq2 != "-1":

                # Only when both HLA alleles can get the corresponding HLA sequence information.

                l_temp = []

                for i in range(0, len(Seq1)):
                    l_temp.append(Seq1[i])
                    l_temp.append(Seq2[i])

                # Reversing
                __return__ = '\t'.join(l_temp)

            else:
                __return__ = '\t'.join(["0" for z in range(0, 2 * _seq_length)])


            return __return__


        else:

            # As each Strings.

            t_Seq1 = ""
            t_Seq2 = ""

            if Seq1 != "-1" and Seq2 != "-1":
                t_Seq1 = Seq1
                t_Seq2 = Seq2
            else:
                t_Seq1 = ''.join(["0" for z in range(0, _seq_length)])
                t_Seq2 = ''.join(["0" for z in range(0, _seq_length)])

            return str([[_HLA_allele1, t_Seq1], [_HLA_allele2, t_Seq2]])




if __name__ == '__main__':

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

    parser.add_argument("-chped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-dict", help="\nHLA dictonary file name(ex. 'HLA_DICTIONARY_AA.txt')\n\n", required=True)
    parser.add_argument("-type", help="\nAA(for Amino Acid) or SNP(for SNPs)\n\n", choices=["AA", "SNPS"], required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as original version.\n\n",
                        action='store_true')
    parser.add_argument("--asLump", help="\n(for Testing) Not zipped result to check strings.\n\n",
                        action='store_true')




    ##### <for Test> #####

    ## (2019. 01. 06.) Introducinig compatibility to work with old version of dictionary

    # # AA
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference_old/HLA_DICTIONARY_AA.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/HAPMAP_CEU.old.enCODED",
    #                           "--previous-version"])

    # # SNPS
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference_old/HLA_DICTIONARY_SNPS.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/HAPMAP_CEU.old.enCODED",
    #                           "--previous-version"])


    ## (2019. 01. 08.) Generalized 4-field

    # # AA
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HAPMAP_CEU_HLA.4field"])

    # # SNPS
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HAPMAP_CEU_HLA.4field"])

    # # AA
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HAPMAP_CEU_HLA.4field.lumped",
    #                           "--asLump"])

    # # SNPS
    # args = parser.parse_args(["-chped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HAPMAP_CEU_HLA.4field.lumped",
    #                           "--asLump"])


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)


    HLAtoSequences(args.chped, args.dict, args.type, args.o, __previous_version=args.previous_version, __asLump=args.asLump)
