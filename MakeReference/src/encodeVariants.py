# -*- coding: utf-8 -*-

import os, re
import argparse, textwrap


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def encodeVariants(_ped, _map, _out, __asSmallLetter=True, __addDummyMarker=False):


    ### Intermediate path.
    _out = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(_out)):
        INTERMEDIATE_PATH = os.path.dirname(_out)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        INTERMEDIATE_PATH = "./"


    # Processing line by line with python Generators to save memory

    ########## < Control Flags > ##########

    _1_ALLELE_OVERLAPPING = 1
    _2_MAKING_NEW_PEDFILE = 1
    _3_MAKING_NEW_MAPFILE = 1
    _4_MAKING_ALLELELIST = 1


    if _1_ALLELE_OVERLAPPING:

        ########## < [1] Allele overlapping > ##########

        # Acquiring column number
        n_loci = 0
        n_row = 0

        with open(_ped, 'r') as f:

            for line in f:

                t_line = re.split(r'\s+', line.rstrip('\n'))
                # print(t_line[6:])

                if n_row == 0:

                    genomic_info = t_line[6:]
                    n_loci = int(len(genomic_info)/2) # == len(l_factors)

                    # Initializing the list containing factors which appear in each locus.
                    l_factors = [[] for i in range(0, n_loci)] # Initialization


                for i in range(0, n_loci):

                    idx1 = 2*i + 6 # index for `t_line`
                    idx2 = idx1 + 1

                    ##### Allele overlapping
                    if t_line[idx1] != "0" and (t_line[idx1] not in l_factors[i]):
                        l_factors[i].append(t_line[idx1])
                    if t_line[idx2] != "0" and (t_line[idx2] not in l_factors[i]):
                        l_factors[i].append(t_line[idx2])


                n_row += 1
                # if n_row > 5 : break


        # Sorting elements of each lists.
        for i in range(0, len(l_factors)):
            l_factors[i].sort()

        ### --- `l_factors` done.




    if _2_MAKING_NEW_PEDFILE:

        ########## < [2] Making new .ped file > ##########

        with open(_out + ".ped", 'w') as f_NewPed:
            f_NewPed.writelines(MakeNewPed(_ped, l_factors, __asSmallLetter, __addDummyMarker))



    if _3_MAKING_NEW_MAPFILE:

        ########## < [3] Making new .map file > ##########

        with open(_out + ".map", 'w') as f_NewMap:
            f_NewMap.writelines(MakeNewMap(_map, l_factors, __addDummyMarker))



    if _4_MAKING_ALLELELIST:

        ########## < [4] Making *.allelelist file > ##########

        with open(_out + ".factors", 'w') as f_allelelist:
            f_allelelist.writelines(MakeAlleleList(_map, l_factors))




def divideToBinaryMarkers(_SNP1, _SNP2, _factors, __asSmallLetter=True):

    _present_ = "p" if __asSmallLetter else "P"
    _absent_ = "a" if __asSmallLetter else "A"

    Seq = []

    if len(_factors) > 2:

        for j in range(0, len(_factors)):

            if _SNP1 == "0" or _SNP2 == "0":
                Seq.append("0"); Seq.append("0")

            else:

                if _factors[j] == _SNP1:
                    Seq.append(_present_)
                else:
                    Seq.append(_absent_)

                if _factors[j] == _SNP2:
                    Seq.append(_present_)
                else:
                    Seq.append(_absent_)

        if len(_factors) > 3:

            j_end = 1 if len(_factors) == 4 else len(_factors)

            for j in range(0, j_end):

                for k in range(j + 1, len(_factors)):

                    if _SNP1 == "0" or _SNP2 == "0":
                        Seq.append("0"); Seq.append("0")

                    else:
                        if _factors[j] == _SNP1 or _factors[k] == _SNP1:
                            Seq.append(_present_)
                        else:
                            Seq.append(_absent_)

                        if _factors[j] == _SNP2 or _factors[k] == _SNP2:
                            Seq.append(_present_)
                        else:
                            Seq.append(_absent_)

            if len(_factors) > 5:

                j_end = 1 if len(_factors) == 6 else len(_factors)

                for j in range(0, j_end):
                    for k in range(j + 1, len(_factors)):
                        for l in range(k + 1, len(_factors)):

                            if _SNP1 == "0" or _SNP2 == "0":
                                Seq.append("0"); Seq.append("0")

                            else:
                                if _factors[j] == _SNP1 or _factors[k] == _SNP1 or _factors[l] == _SNP1:
                                    Seq.append(_present_)
                                else:
                                    Seq.append(_absent_)

                                if _factors[j] == _SNP2 or _factors[k] == _SNP2 or _factors[l] == _SNP2:
                                    Seq.append(_present_)
                                else:
                                    Seq.append(_absent_)



    else:
        # Most of Cases have length less than equal 2, they will fall into this if-else block.
        if _SNP1 == "0" or _SNP2 == "0":
            Seq.append("0"); Seq.append("0")
        else:
            Seq.append(_SNP1); Seq.append(_SNP2)

    return '\t'.join(Seq)



def MakeNewPed(_p_ped, _l_factors, __asSmallLetter=True, __addDummyMarker=False):

    count = 0

    with open(_p_ped, 'r') as f:
        for line in f:
            t_line = re.split(r'\s+', line.rstrip('\n'))

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_info__ = '\t'.join([
                divideToBinaryMarkers(t_line[2 * i + 6], t_line[2 * i + 7], _l_factors[i], __asSmallLetter) for i in range(0, len(_l_factors))
            ])

            __return__ = '\t'.join([__ped_info__, __genomic_info__])


            if __addDummyMarker:
                # add Dummy Markers.
                dummy_markers = '\t'.join(['d', 'D'] if bool(count % 2) else ['D', 'd'])
                __return__ = '\t'.join([__return__, dummy_markers])


            yield __return__ + "\n"




def MakeNewMap(_p_map, _l_factors, __addDummyMarker=False):

    count = 0

    with open(_p_map, 'r') as f_map:

        for l in f_map:

            idx = count # index for `l_factors`

            if len(_l_factors[idx]) > 2:

                t_line = re.split(r'\s+', l.rstrip('\n'))

                for j in range(0, len(_l_factors[idx])):
                    yield '\t'.join([t_line[0], t_line[1] + '_' + _l_factors[idx][j], t_line[2], t_line[3]]) + "\n"


                if len(_l_factors[idx]) > 3:

                    j_end = 1 if len(_l_factors[idx]) == 4 else len(_l_factors[idx])

                    for j in range(0, j_end):
                        for k in range(j+1, len(_l_factors[idx])):
                            yield '\t'.join([t_line[0], t_line[1]+ '_' + _l_factors[idx][j]+_l_factors[idx][k], t_line[2], t_line[3]]) + "\n"


                    if len(_l_factors[idx]) > 5:

                        j_end = 1 if len(_l_factors[idx]) == 6 else len(_l_factors[idx])

                        for j in range(0, j_end):
                            for k in range(j+1, len(_l_factors[idx])):
                                for l in range(k+1, len(_l_factors[idx])):
                                    yield '\t'.join([t_line[0], t_line[1]+ '_' + _l_factors[idx][j]+_l_factors[idx][k]+_l_factors[idx][l], t_line[2], t_line[3]]) + "\n"

            else:
                yield l # "\n" is included.

            count += 1


        if __addDummyMarker:
            # Adding Dummy Marker
            yield '\t'.join(["6", "dummy_marker", "0", "33999999"]) + "\n"



def MakeAlleleList(_p_map, _l_factors):

    count = 0

    with open(_p_map, 'r') as f_map:

        for l in f_map:

            idx = count

            t_line = re.split(r'\s+', l.rstrip('\n'))

            locus_label = t_line[1]
            alleleset = _l_factors[idx]

            yield '\t'.join([locus_label, str(alleleset)]) + "\n"

            count += 1




if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        Original Author : Sherman Jia, 2012
     
        encodeVariants.py
     
        - This script generates PLINK binary markers which encodes multi-alleles(factors) at 
            a locus (Amino acids and SNPs) 

    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nThe *.ped file which is generated by 'HLAtoSequences.py'.\n\n", required=True)
    parser.add_argument("-map", help="\nThe *.map file which is generated by 'HLAtoSequences.py'.\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)


    ##### <for Test> #####

    # (2018. 7. 13.)
    # # AA
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.enCODED"])

    # # SNPS
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.enCODED"])

    # (2019. 1. 6.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/_Case_HAPMAP_CEU.old.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/HLA_DICTIONARY_AA.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_2_encodeVariants/_Case_HAPMAP_CEU.AA.enCODED"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/_Case_HAPMAP_CEU.old.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/HLA_DICTIONARY_SNPS.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_2_encodeVariants/_Case_HAPMAP_CEU.SNPS.enCODED"])

    # (2019. 1. 9.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HAPMAP_CEU_HLA.4field.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HLA_DICTIONARY_AA.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_2_encodeVariants/HAPMAP_CEU_HLA.4field.AA.enCODED"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HAPMAP_CEU_HLA.4field.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_1_HLAtoSequences/HLA_DICTIONARY_SNPS.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190109/_2_encodeVariants/HAPMAP_CEU_HLA.4field.SNPS.enCODED"])


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    encodeVariants(args.ped, args.map, args.o)