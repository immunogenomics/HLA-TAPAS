# -*- coding: utf-8 -*-

# (2017/11/27) recoded by Wanson Choi

###########################################################################################
#
# Sherman Jia, 2012
# encodeHLA.pl
#
# This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
# The input ped file should contain: FID,IID,pID,mID,SEX,PHENO,
#                                    2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1
#
###########################################################################################


# Breifly, This module will generate a new .map and .ped file which have information of HLA alleles marked as 'A' or 'P'.

import os, re
import pandas as pd
import argparse, textwrap
from collections import OrderedDict


def encodeHLA(_INPUT_PED, _OUTPUT, _hg = "19"):


    ########## <Core Variables> ##########

    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

    # isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}

    # (2018. 3. 22.) Modify to use "HLA_INTEGRATED_POSITIONS_hg{18,19,38}.txt"
    # (2018. 7. 16.) Statically prepared beforehand. Originally, it is supposed to utilize "HLA_EXON_POSITIONS_SNPS_hg{18,19,38}.txt" file.
    genepos_hg ={"18" : {"A": 30018226, "C": 31344505, "B": 31429628, "DRB1": 32654525, "DQA1": 32713161, "DQB1": 32735219, "DPA1": 33140324, "DPB1": 33151681},
                 "19" : {"A": 29910247, "C": 31236526, "B": 31321643, "DRB1": 32546547, "DQA1": 32605183, "DQB1": 32627241, "DPA1": 33032346, "DPB1": 33043703},
                 "38" : {"A": 29942470, "C": 31268749, "B": 31353866, "DRB1": 32578770, "DQA1": 32637406, "DQB1": 32659464, "DPA1": 33064569, "DPB1": 33075926}}

    # Module Name for stdout
    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME + "Init.")


    # 주어진 ped파일에 HLA column 별로 나타나는 allele name들의 집합
    ALLELE_TABLES = OrderedDict()
    ALLELE_TABLES_1field = OrderedDict()

    # map파일의 Lable로 준비시킬 항목들 집합
    ALL_ALLELES = [] # 결국 ALL_ALLELES == map_LABELS 라고 생각해도 괜춘.
    dict_ALL_ALLELES = {} # 추후에 search할때 빨리 하려고


    df_OUTPUT_map = pd.DataFrame()
    df_OUTPUT_ped = pd.DataFrame()




    ########## < Control Flags > ##########


    LOADING_PEDFILE = 1
    MAKING_ALLELE_TABLE = 1
    MAKING_OUTPUT_MAP = 1
    MAKING_OUTPUT_PED = 1





    if LOADING_PEDFILE:

        ########## < 1. Loading Input PED file > ##########

        print(std_MAIN_PROCESS_NAME + "[1] Loading Input PED file.\n")

        INPUT_PED = pd.read_table(_INPUT_PED, sep='\t', header=None, dtype=str,
                                  names = ['Fam_ID', 'Sample_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phe'] + [''.join([HLA_names[i], '_', str(j)]) for i in range(0, len(HLA_names)) for j in range(1,3)],
                                  index_col=[0, 1, 2, 3, 4, 5]
                                  )

        print(INPUT_PED.head())



    if MAKING_ALLELE_TABLE:

        ########## < 2. Making Allele Table > ##########

        print(std_MAIN_PROCESS_NAME + "[2] Making Allele Table.\n")

        """
        for i in range(0, len(INPUT_PED.index)):

            line = tuple(INPUT_PED.iloc[i, :])
            # 0:6(6컬럼) => .ped sample header
            # 6: => HLA allele pair

            alleles["HLA_A_"+line[6]] = genepos["HLA_A"];       alleles["HLA_A_"+line[7]] = genepos["HLA_A"];
            alleles["HLA_A_"+line[6][0:2]] = genepos["HLA_A"];  alleles["HLA_A_"+line[7][0:2]] = genepos["HLA_A"];

            alleles["HLA_B_"+line[8]] = genepos["HLA_B"];       alleles["HLA_B_"+line[9]] = genepos["HLA_B"];
            alleles["HLA_B_"+line[8][0:2]] = genepos["HLA_B"];  alleles["HLA_B_"+line[9][0:2]] = genepos["HLA_B"];

            alleles["HLA_C_"+line[10]] = genepos["HLA_C"];       alleles["HLA_C_"+line[11]] = genepos["HLA_C"];
            alleles["HLA_C_"+line[10][0:2]] = genepos["HLA_C"];  alleles["HLA_C_"+line[11][0:2]] = genepos["HLA_C"];

            alleles["HLA_DPA1_"+line[12]] = genepos["HLA_DPA1"];       alleles["HLA_DPA1_"+line[13]] = genepos["HLA_DPA1"];
            alleles["HLA_DPA1_"+line[12][0:2]] = genepos["HLA_DPA1"];  alleles["HLA_DPA1_"+line[13][0:2]] = genepos["HLA_DPA1"];

            alleles["HLA_DPB1_"+line[14]] = genepos["HLA_DPB1"];       alleles["HLA_DPB1_"+line[15]] = genepos["HLA_DPB1"];
            alleles["HLA_DPB1_"+line[14][0:2]] = genepos["HLA_DPB1"];  alleles["HLA_DPB1_"+line[15][0:2]] = genepos["HLA_DPB1"];

            alleles["HLA_DQA1_"+line[16]] = genepos["HLA_DQA1"];       alleles["HLA_DQA1_"+line[17]] = genepos["HLA_DQA1"];
            alleles["HLA_DQA1_"+line[16][0:2]] = genepos["HLA_DQA1"];  alleles["HLA_DQA1_"+line[17][0:2]] = genepos["HLA_DQA1"];

            alleles["HLA_DQB1_"+line[18]] = genepos["HLA_DQB1"];       alleles["HLA_DQB1_"+line[19]] = genepos["HLA_DQB1"];
            alleles["HLA_DQB1_"+line[18][0:2]] = genepos["HLA_DQB1"];  alleles["HLA_DQB1_"+line[19][0:2]] = genepos["HLA_DQB1"];

            alleles["HLA_DRB1_"+line[20]] = genepos["HLA_DRB1"];       alleles["HLA_DRB1_"+line[21]] = genepos["HLA_DRB1"];
            alleles["HLA_DRB1_"+line[20][0:2]] = genepos["HLA_DRB1"];  alleles["HLA_DRB1_"+line[21][0:2]] = genepos["HLA_DRB1"];


        for k,v in alleles.items():
            print("key : {0} / value : {1}".format(k, v))
            
        """

        # In the past, `ALLELE_TABLES` was created based on dictionary data structure, but now it is created by "apply()" function with "set()" function.

        for i in range(0, len(HLA_names)):
            # for i in range(0, 1):

            temp = INPUT_PED.filter(regex=HLA_names[i]+'_\d', axis=1).apply(set, axis=0).apply(lambda x : x.difference({0, "0"}))

            # print(temp)

            # Column 1 and 2
            set_al1 = temp.iat[0] # ex)     A_1 := {A*32:01:01, A*23:01:01, A*26:01:01, A*02:01:0...
            set_al2 = temp.iat[1] # ex)     A_2 := {A*24:02:01:01, A*01:01:01:01, A*02:05:01, A*3...

            # sr_Unioned_Set = pd.Series(list(set_al1.union(set_al2))).sort_values() # sorting

            l_Unioned_Set = list(set_al1.union(set_al2))
            l_Unioned_Set.sort() # sorting

            # l_Unioned_Set = pd.Series(l_Unioned_Set)
            print(l_Unioned_Set)

            ALLELE_TABLES[HLA_names[i]] = l_Unioned_Set


            ##### Dealing with 1-field #####

            if len(l_Unioned_Set) > 0:

                ### The case where the union of set of each two column has at least 1 element.

                sr_temp_1field = pd.Series(l_Unioned_Set).apply(lambda x : re.match(pattern='\*'.join([HLA_names[i], '\d{2,3}']), string=x).group()).unique() # Going through "unique()" function.

                # print("\nsr_temp_1field\n")
                # print(sr_temp_1field)

                ALLELE_TABLES_1field[HLA_names[i]] = sr_temp_1field.tolist()

            else:
                ### 집합 원소의 개수가 0 일때,(아예 input_ped에서 allele_name이 '0'으로 주어져서 없는 경우)
                ALLELE_TABLES_1field[HLA_names[i]] = l_Unioned_Set


        # for k, v in ALLELE_TABLES.items():
        #
        #     print("\n===============\n")
        #     print("{0} : \n{1}".format(k, v))



    if MAKING_OUTPUT_MAP:

        ########## < 3. Making OUTPUT .map file > ##########

        print(std_MAIN_PROCESS_NAME + "[3] Making OUTPUT .map file.\n")

        """        
        to_df_OUTPUT_map = []


        for name in HLA_names:
            for k in sorted_keys:
                temp = k.split('_')
                al = temp[2]

                if (temp[1] == name) and (al != "NA") and (al != "") and (al != "0") and (al != "0 0"):
                    # to_df_OUTPUT_map += [ ["6", k, "0", genepos['_'.join([temp[0], temp[1]])]] ]
                    to_df_OUTPUT_map.extend([("6", k, "0", genepos['_'.join([temp[0], temp[1]])])])

        df_OUTPUT_map = pd.DataFrame(to_df_OUTPUT_map)
        df_OUTPUT_map.to_csv(_OUTPUT + '.HLA.map', sep='\t', header=False, index=False)
                
        """


        """ 
        map_Label 만들기 
        
        (1) map_LABELS
        (2) map_CHR
        (3) map_GENETIC_DISTANCE
        (4) map_POS
        
        """

        ##### Making Label for *.map file. #####

        # ALL_ALLELES = pd.concat([ALLELE_TABLES[HLA_names[i]].append(ALLELE_TABLES_1field[HLA_names[i]]) for i in range(0, len(HLA_names))]).sort_values()

        ALL_ALLELES = ALLELE_TABLES[HLA_names[0]]
        ALL_ALLELES.extend(ALLELE_TABLES_1field[HLA_names[0]])
        # ALL_ALLELES = [ALLELE_TABLES[HLA_names[i]].append(ALLELE_TABLES_1field[HLA_names[i]]) for i in range(0, len(HLA_names))]

        for i in range(1, len(HLA_names)):
            ALL_ALLELES.extend(ALLELE_TABLES[HLA_names[i]])
            ALL_ALLELES.extend(ALLELE_TABLES_1field[HLA_names[i]])

        ALL_ALLELES.sort()

        print("\nALL_ALLELES\n")
        print(ALL_ALLELES)


        ### HLA_index(Mining HLA gene name)
        sr_HLA = pd.Series(ALL_ALLELES).apply(lambda x : re.search(pattern='\w+\*', string=x).group().rstrip('*'))
        # print(sr_HLA)

        ### map_LABELS & map_POS ###
        map_LABELS = pd.Series(['HLA_' + ALL_ALLELES[i] for i in range(0, len(ALL_ALLELES))])
        # print("\n`map_LABELS`\n")
        # print(map_LABELS)

        map_POS = [genepos_hg[_hg][sr_HLA.iat[i]] for i in range(0, len(sr_HLA))]
        # print(map_POS)


        # 추후 ped파일에서 search할때 편하게 HLA_name별로 나눠놓기
        dict_ALL_ALLELES = {HLA_names[i] : [ALL_ALLELES[j] for j in range(0, len(ALL_ALLELES)) if (HLA_names[i]+'*' in ALL_ALLELES[j])] for i in range(0, len(HLA_names))}

        # print("\nsegmented `ALL_ALLELES`\n")
        # for k,v in dict_ALL_ALLELES.items():
        #     print("\n============\n")
        #     print("{0} : \n{1}".format(k, v))


        ##### map_Label을 제외한 나머지 map파일 항목 만들기 #####

        map_CHR = ['6' for i in range(0, len(map_LABELS))]
        map_GENETIC_DISTANCE = ['0' for i in range(0, len(map_LABELS))]


        df_OUTPUT_map = pd.DataFrame.from_dict({"Chr" : map_CHR, "Name" : map_LABELS.tolist(), "GD" : map_GENETIC_DISTANCE, "POS" : map_POS}).loc[:, ["Chr", "Name", "GD", "POS"]]
        print(std_MAIN_PROCESS_NAME + "Output .map file.\n")
        print(df_OUTPUT_map.head(50))

        df_OUTPUT_map.to_csv('.'.join([_OUTPUT,'HLA.map']), sep='\t', header=False, index=False)



    if MAKING_OUTPUT_PED:

        ########## < 4. Making OUTPUT.ped file > ##########

        print(std_MAIN_PROCESS_NAME + "[4] Making .ped file.\n")


        """
        
                to_df_OUTPUT_ped = []
        
                for i in range(0, len(INPUT_PED.index)):
        
                    line = tuple(INPUT_PED.iloc[i, :])
        
                    to_df_OUTPUT_ped.extend([
                        line[0:6] +
                        PrintGenotypes("A", line[6], line[7], sorted_keys) +
                        PrintGenotypes("C", line[10], line[11], sorted_keys) +
                        PrintGenotypes("B", line[8], line[9], sorted_keys) +
                        PrintGenotypes("DRB1", line[20], line[21], sorted_keys) +
                        PrintGenotypes("DQA1", line[16], line[17], sorted_keys) +
                        PrintGenotypes("DQB1", line[18], line[19], sorted_keys) +
                        PrintGenotypes("DPA1", line[12], line[13], sorted_keys) +
                        PrintGenotypes("DPB1", line[14], line[15], sorted_keys)
                    ])
        
                # 실제로 HLA allele들이 genome상에서 저 position순으로 존재하기 때문에 저렇게 구성시키는게 맞다고함
                # 가지고 시작하는 HAPMAP_CEU_HLA.ped파일의 HLA allele순서 구성은 정말 우리들끼리 임의대로 그렇게 구성해놓은 것일 뿐임.
                # from 승호쌤's comment
                
                df_OUTPUT_ped = pd.DataFrame(to_df_OUTPUT_ped)
                df_OUTPUT_ped.to_csv(_OUTPUT + '.HLA.ped', sep='\t', header=False, index=False)
                    
                """


        to_df_OUTPUT_ped = []

        # for i in range(0, 5):
        for i in range(0, INPUT_PED.shape[0]):

            # print("\n================\n")

            line_INPUT_PED = tuple(INPUT_PED.iloc[i, :])
            # print(line_INPUT_PED)


            t_line_OUTPUT_PED = [PrintGenotypes3(line_INPUT_PED[2*j], line_INPUT_PED[2*j+1], dict_ALL_ALLELES[HLA_names[j]]) for j in range(0, len(HLA_names))]
            # print(t_line_OUTPUT_PED)
            # print(pd.Series(dict_ALL_ALLELES["A"]))

            # Flattening
            line_OUTPUT_PED = [item for eachlist in t_line_OUTPUT_PED for item in eachlist]

            # print("\nFlattened t_line_OUTPUT_PED is \n")
            # print(line_OUTPUT_PED)

            to_df_OUTPUT_ped.append(line_OUTPUT_PED)


        df_OUTPUT_ped = pd.DataFrame(to_df_OUTPUT_ped)
        df_OUTPUT_ped.index = INPUT_PED.index

        print(df_OUTPUT_ped.head())


        df_OUTPUT_ped.to_csv('.'.join([_OUTPUT, 'HLA.ped']), sep='\t', header=False, index=True)



    return 0


def PrintGenotypes3(_allele1, _allele2, _seg_ALL_ALLELES):

    l_output = []

    # print("\nAlleles : {0} and {1}".format(_allele1, _allele2))
    # print("\ndict_ALL_ALLELES: \n{0}".format(_seg_ALL_ALLELES))

    if len(_seg_ALL_ALLELES) > 0:

        if ((_allele1 != "0") and (_allele2 != "0")):  # The condition for checking integer value 0 won't be included here because .ped file was read with "dtype=str" option.

            t_sr1 = pd.Series(_seg_ALL_ALLELES, index=pd.Index(_seg_ALL_ALLELES)).apply(lambda x: "P" if (x in _allele1) else "A")
            # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr1))
            t_sr2 = pd.Series(_seg_ALL_ALLELES, index=pd.Index(_seg_ALL_ALLELES)).apply(lambda x: "P" if (x in _allele2) else "A")
            # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr2))

            for i in range(0, len(t_sr1)):
                l_output.append(t_sr1.iat[i])
                l_output.append(t_sr2.iat[i])

            # print(l_output)
            return l_output


        else:

            # Comments... at most below part.
            # print("At least one of allele is 0")

            for i in range(0, len(_seg_ALL_ALLELES)):
                l_output.append("0")
                l_output.append("0")

            return l_output

    else:
        # In cases such as "DPA1" or "DPB1 where any alleles don't appear, Just skip.
        # Then just return NULL
        return l_output



if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        Sherman Jia, 2012
        encodeHLA.py

        This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
        The input ped file should contain: FID,IID,pID,mID,SEX,PHENO,
                                            2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1

    ###########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")



    ##### <for Test> #####

    # (2018.2.28 모듈 테스트)
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./TEST_0228", "--hg", "19"])

    # args = parser.parse_args(["./HAPMAP_CEU_HLA.ped", "./TEST_0228", "--hg", "19"])
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./TEST_0305", "-hg", "19"])

    # (2018. 7. 16.)

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-hg", "19",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370"])



    ##### <for Publication> #####

    args = parser.parse_args()


    print(args)

    # Implementing Main Function
    encodeHLA(args.ped, args.o, args.hg)




#========== (deprecated) ==========

# def PrintGenotypes(_locus, _allele1, _allele2, _sorted_keys):
#
#     G1 = "0"
#     G2 = "0"
#
#     l_output = []
#
#     for l in _sorted_keys:
#
#         temp = l.split('_')
#         temp_allele = temp[2]
#         temp_locus = temp[1]
#
#         # (2018. 3. 7)
#         # ex) 예전꺼 기준으로 _sorted_keys에는 'HLA_A_01', 'HLA_A_0101', 'HLA_A_02', 'HLA_A_0201', 'HLA_A_0203', 'HLA_A_0205', ... 들이 들어있음.
#         # 저렇게 '_'로 split하게 되면 temp_allele에는 '01', '0101'등이 들어가있는거고, temp_locus에는 'A', 'C', 'B' 등이 들어가게 되는 거임.
#         # 지금은 4-field name system으로 일반화된 이름들이 ALLELE_TABLES dictionary에 담기게 됨.
#
#
#         if (temp_allele != "NA" and temp_allele != "" and temp_allele != "0" and temp_allele != "0 0") and temp_locus == _locus:
#
#             if ((_allele1 != "NA" and _allele1 != "" and _allele1 != "0") and
#                 (_allele2 != "NA" and _allele2 != "" and _allele2 != "0")):
#
#                 # allele1 에 대한 genotype 값 구하기("P" or "A")
#                 if _allele1 == "NA" or _allele1 == "" or _allele1 == "0":
#                     G1 = "0"
#                 else:
#                     if _allele1 == temp_allele or _allele1[0:2] == temp_allele:
#                         G1 = "P"
#                     elif (len(_allele1) == 2 and len(temp_allele) == 4 and temp_allele[0:2] == _allele1):
#                         G1 = "0"
#                     else:
#                         G1 = "A"
#
#
#                 # allele2에 대한 genotype 값 구하기
#                 if _allele2 == "NA" or _allele2 == "" or _allele2 == "0":
#                     G2 = "0"
#                 else:
#                     if _allele2 == temp_allele or _allele2[0:2] == temp_allele:
#                         G2 = "P"
#                     elif (len(_allele2) == 2 and len(temp_allele) == 4 and temp_allele[0:2] == _allele2):
#                         G2 = "0"
#                     else:
#                         G2 = "A"
#
#                 if G1 == "0" or G2 == "0":
#                     # l_output += ["0", "0"]
#                     l_output.extend(("0", "0"))
#                 else:
#                     # l_output += [G1, G2]
#                     l_output.extend((G1, G2))
#
#
#             else:
#                 # l_output += ["0", "0"]
#                 l_output.extend(("0", "0"))
#
#     return tuple(l_output)


# def PrintGenotypes2(_locus, _allele1, _allele2, _seg_ALL_ALLELES):
#
#
#     l_output = []
#
#     print("\nAlleles : {0} and {1}".format(_allele1, _allele2))
#
#     """
#     이 함수의 본질은 결국 map파일의 Label의 항목들(ALL_ALLELES, 더 정확히는 seg_ALL_ALLELES)에 대해서 한번씩 for문을 돈다는 거임.
#     """
#
#     # 어차피 주어진 allele들 중 둘 중 하나라도 삐꾸면 어차피 ["0", "0"] 만 오질라게 append시켜야함.
#
#     if len(_seg_ALL_ALLELES[_locus]) > 0:
#
#         if ((_allele1 != "0") and
#             (_allele2 != "0")):     # input ped파일 읽어들일때 dtype=str했으니 숫자 0이랑 비교하는 조건식은 빼겠음.
#
#             t_sr1 = pd.Series(_seg_ALL_ALLELES[_locus], index=pd.Index(_seg_ALL_ALLELES[_locus])).apply(lambda x : "P" if (x in _allele1) else "A")
#             # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr1))
#             t_sr2 = pd.Series(_seg_ALL_ALLELES[_locus], index=pd.Index(_seg_ALL_ALLELES[_locus])).apply(lambda x : "P" if (x in _allele2) else "A")
#             # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr2))
#
#             for i in range(0, len(t_sr1)):
#                 l_output.append(t_sr1.iat[i])
#                 l_output.append(t_sr2.iat[i])
#
#             # print(l_output)
#             return l_output
#
#
#         else:
#
#             # (2018.4.2) 여기만 다시 확인해볼것. 원래 소스에서 DPA1, DPB1같이 원래 값 0ㅇ인애들은 어떻게 escape시켰는지 다시 확인좀.
#
#             """
#             사실 생각해보면 그 문제인거같음.
#             input으로 주어진 ped파일에, HLA 컬럼별로 한 쌍의 allele이 주어졌을 때 여기서 하나라도 0인 애들은 사실 map파일 Label집합에도 포함되지 않음. 구체적으로는 MAKING_ALLELE_TABLE 코드블럭에서
#
#                 temp = INPUT_PED.filter(regex=HLA_names[i]+'_\d', axis=1).apply(set, axis=0).apply(lambda x : x.difference({0, "0"}))
#
#             이 코드 부분에서 뒤쪽에 0인애들은 set_difference시켜버렸으니 확실히 map_Label 집합에는 DPA1, DPB1같은 애들이 포함되지는 않을 것임. 결과적으로 seg_ALL_ALLELE["DPA1"] == [] 임을 보장해 준다는 것.
#
#             정리하자면, 현재 if-else블럭 밖에 len(seg_ALL_ALLELES[HLA_name[i]]) > 0 인 경우로 한번만 더 감싸면 DPA1, DPB1같은 애들을 아예 고려하지 않는 기존의 encodeHLA.pl스럽게 할 수 있는거 아니냐는 거지.
#
#             """
#
#
#             print("At least one of allele is 0")
#
#             for i in range(0, len(_seg_ALL_ALLELES[_locus])):
#                 l_output.append("0")
#                 l_output.append("0")
#
#             return l_output
#
#     else:
#         # DPA1, DPB1처럼 아예 allele이 나타나지 않은 경우는 그냥 스킵(어차피 map파일의 label집합상에 안나타났을테니
#         # 이럴때는 그냥 NULL list 리턴(어차피 flattening하는 과정에서 없애줌.
#         return l_output


# (2018.4.2) 여기만 다시 확인해볼것. 원래 소스에서 DPA1, DPB1같이 원래 값 0ㅇ인애들은 어떻게 escape시켰는지 다시 확인좀.

"""
사실 생각해보면 그 문제인거같음.
input으로 주어진 ped파일에, HLA 컬럼별로 한 쌍의 allele이 주어졌을 때 여기서 하나라도 0인 애들은 사실 map파일 Label집합에도 포함되지 않음. 구체적으로는 MAKING_ALLELE_TABLE 코드블럭에서

    temp = INPUT_PED.filter(regex=HLA_names[i]+'_\d', axis=1).apply(set, axis=0).apply(lambda x : x.difference({0, "0"}))

이 코드 부분에서 뒤쪽에 0인애들은 set_difference시켜버렸으니 확실히 map_Label 집합에는 DPA1, DPB1같은 애들이 포함되지는 않을 것임. 결과적으로 seg_ALL_ALLELE["DPA1"] == [] 임을 보장해 준다는 것.

정리하자면, 현재 if-else블럭 밖에 len(seg_ALL_ALLELES[HLA_name[i]]) > 0 인 경우로 한번만 더 감싸면 DPA1, DPB1같은 애들을 아예 고려하지 않는 기존의 encodeHLA.pl스럽게 할 수 있는거 아니냐는 거지.

"""
