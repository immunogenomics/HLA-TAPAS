# -*- coding: utf-8 -*-


import os, sys, re
import argparse, textwrap
import pandas as pd


def ClassifyGroups(_p_allelelist, _p_Ggroup, _p_Pgroup, _p_dict_AA, _p_dict_SNPS, _out, _imgt = ""):

    """

    This function takes "hla_nom_g.txt" or "hla_nom_p.txt" file and converts them to
    `INTEGRATED_ALLELELIST` file(".iat" file.)

    The output will be majorly used in "COATING_PED.py" module.

    Sometimes there is no matched group for an allele especially in case of P-group.
    In this case the matched value will be filled with -1 which represents "Not Found".

    """


    ########## < Core Variables > ##########

    HLA_names  = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME+"Init.")

    # set_al_Allelelist = None
    HLA_DICT_AA = None
    HLA_DICT_SNPS = None


    ##### < Control Flags > #####

    _1_LOADING_ALLELELIST = 1
    _2_READING_HLA_NOM_G = 1
    _3_READING_HLA_NOM_P = 1
    _4_LOADING_HLA_DICTIONARY = 1
    _5_MERGE_TABLES = 1
    _6_MAIN_MATCHING = 1



    if _1_LOADING_ALLELELIST:

        ########## < 1. Reading "Allelelist.txt" file > ##########

        header_lines = []
        major_lines = []

        f = open(_p_allelelist, 'r')
        temp_string = f.readline()

        if temp_string.startswith('#'):

            # I don't know why but "Allelelist.txt" file's header and separter changed between below-imgt3320 and above-equal-imgt3320.

            header_lines.append(temp_string)

            while(temp_string.startswith('#')):

                temp_string = f.readline()
                header_lines.append(temp_string)

                # Not intended, but anyway the line which to be column index will be includes as last element in `header` list('AlleleID,Allele\n').

            major_lines = f.readlines()

        else:

            major_lines.append(temp_string)
            major_lines.extend(f.readlines())

        print("\nheader lines are:\n{0}".format(header_lines))
        print("\nMajor lines are:\n{0}".format(major_lines))

        # Check delimeter
        p = re.compile(r'\s+')
        p_delimeter = re.compile(r'\s+') if bool(p.search(major_lines[0].rstrip('\n'))) else re.compile(',')

        df_Allelelist = pd.DataFrame([p_delimeter.split(item.rstrip('\n')) for item in major_lines],
                                     columns = pd.Index(["AlleleID", "Allele"]))

        print(std_MAIN_PROCESS_NAME+"Loaded Allelelist file.\n")
        print(df_Allelelist.head())

        # (2018. 7. 15.) Filtering only major HLA genes in `df_Allelelist` beforehand.

        p_Major_HLA_Genes = re.compile('|'.join([name + "\*" for name in HLA_names]))

        flag1 = df_Allelelist.loc[:, "Allele"].str.match(p_Major_HLA_Genes)

        df_Allelelist = df_Allelelist.loc[flag1]



    if _2_READING_HLA_NOM_G:

        ########## < 2. Reading "hla_nom_g.txt" file > ##########

        f = open(_p_Ggroup, 'r')
        lines = f.readlines()

        idx_start = 0
        while lines[idx_start].startswith('#'): idx_start += 1

        # (2018. 7. 12.) No version checking anymore. In case of imgt370 version, there is no information of version.

        major_lines = lines[idx_start:]

        sr_hla_g = pd.Series(major_lines).apply(lambda x: x.rstrip('\n')).apply(lambda x: x.split(';'))
        df_hla_g = pd.DataFrame(sr_hla_g.tolist(), columns=["HLA", "alleles", "Group"])

        # column trimming
        df_hla_g.loc[:, "HLA"] = df_hla_g.loc[:, "HLA"].apply(lambda x : x.rstrip('*'))
        df_hla_g.loc[:, "alleles"] = df_hla_g.loc[:, "alleles"].apply(lambda x : x.split('/'))

        # setting index with "HLA" column values
        df_hla_g = df_hla_g.set_index("HLA")

        print(std_MAIN_PROCESS_NAME+"Loaded g group file.\n")
        print(df_hla_g.head())

        dict_hla_g = {HLA_names[i] : df_hla_g.loc[HLA_names[i], :] for i in range(0, len(HLA_names))}



    if _3_READING_HLA_NOM_P:

        ########## < 3. Reading "hla_nom_p.txt" file > ##########

        f = open(_p_Pgroup, 'r')
        lines = f.readlines()

        idx_start = 0
        while lines[idx_start].startswith('#'): idx_start += 1

        major_lines = lines[idx_start:]

        sr_hla_p = pd.Series(major_lines).apply(lambda x: x.rstrip('\n')).apply(lambda x: x.split(';'))
        df_hla_p = pd.DataFrame(sr_hla_p.tolist(), columns=["HLA", "alleles", "Group"])

        # column trimming
        df_hla_p.loc[:, "HLA"] = df_hla_p.loc[:, "HLA"].apply(lambda x : x.rstrip('*'))
        df_hla_p.loc[:, "alleles"] = df_hla_p.loc[:, "alleles"].apply(lambda x : x.split('/'))

        # setting index with "HLA" column values
        df_hla_p = df_hla_p.set_index("HLA")

        print(std_MAIN_PROCESS_NAME+"Loaded p group file.\n")
        print(df_hla_p.head())

        dict_hla_p = {HLA_names[i] : df_hla_p.loc[HLA_names[i], :] for i in range(0, len(HLA_names))}



    if _4_LOADING_HLA_DICTIONARY:

        ########## < 4. Loading "HLA_DICTIONARY_{AA,SNPS}.txt" > ##########

        HLA_DICT_AA = pd.read_table(_p_dict_AA, sep='\t', header=None, names=["Allele", "Seqs"])
        HLA_DICT_SNPS = pd.read_table(_p_dict_SNPS, sep='\t', header=None, names=["Allele", "Seqs"])

        print("\nset_al_HLA_DICT_AA\n")
        print(HLA_DICT_AA.head())

        print("\nset_al_HLA_DICT_SNPS\n")
        print(HLA_DICT_SNPS.head())



    if _5_MERGE_TABLES:

        ########## < 5. Merge > ##########

        print(len(df_Allelelist))

        merged_AA = pd.merge(df_Allelelist, HLA_DICT_AA, how='left', on="Allele").fillna("-1")
        merged_SNPS = pd.merge(df_Allelelist, HLA_DICT_SNPS, how='left', on="Allele").fillna("-1")

        print("\nMerged Results 1\n")
        print(merged_AA.head(20))
        # print(len(merged_AA))

        print("\nMerged Results 2\n")
        print(merged_SNPS.head(20))
        # print(len(merged_SNPS))

        SeqisAvailable_AA = merged_AA.loc[:, "Seqs"].apply(lambda x : "P" if x != "-1" else "A")
        SeqisAvailable_SNPS = merged_SNPS.loc[:, "Seqs"].apply(lambda x : "P" if x != "-1" else "A")

        # print("\nLength : {0} and {1}".format(len(SeqisAvailable_AA), len(SeqisAvailable_SNPS)))


        df_Allelelist_SeqInfo = pd.concat([SeqisAvailable_AA, SeqisAvailable_SNPS], axis=1, keys=["AA_Seq", "SNPS_Seq"])

        print(df_Allelelist_SeqInfo.head())



    if _6_MAIN_MATCHING:

        ########## < Matching Job > ##########

        # Main iteration will be conducted on the elements of "alleles" column in `df_Allelelist` dataframe.
        # Before performing it, `df_Allelelist` dataframe should filter alleles of minor HLA genes such as "W*0101".

        df_Allelelist_filtered = pd.DataFrame(df_Allelelist.loc[:, "Allele"].apply(lambda x : x.split('*')).tolist(), columns=["HLA", "Pure_Alleles"], index=pd.Index( df_Allelelist.loc[:, "Allele"]))
        print(df_Allelelist_filtered.head())


        Found_Gs = pd.Series([whichGroup(df_Allelelist_filtered.iat[i, 1], dict_hla_g[df_Allelelist_filtered.iat[i, 0]]) for i in range(0, len(df_Allelelist_filtered))])
        print(Found_Gs)
        # Found_Gs.to_csv("./Found_Gs.txt")

        Found_Ps = pd.Series([whichGroup(df_Allelelist_filtered.iat[i, 1], dict_hla_p[df_Allelelist_filtered.iat[i, 0]]) for i in range(0, len(df_Allelelist_filtered))])
        print(Found_Ps)
        # Found_Ps.to_csv("./Found_Ps.txt")


        print(len(df_Allelelist_SeqInfo))
        print(len(Found_Gs))

        df_OUTPUT = pd.concat([Found_Gs, Found_Ps, df_Allelelist_SeqInfo], axis=1)
        df_OUTPUT.index = df_Allelelist_filtered.index
        df_OUTPUT.columns = ["G_group", "P_group", "AA_Seq", "SNPS_Seq"]

        print(std_MAIN_PROCESS_NAME+"Output of 'i'ntegrated 'a'llele 't'able(.iat).\n")
        print(df_OUTPUT.head(50))

        df_OUTPUT.to_csv('.'.join([_out + ("imgt"+_imgt if bool(_imgt) else ""), "iat"]), sep='\t', header=True, index=True)


        print(std_MAIN_PROCESS_NAME+"Finished.\n")

    return 0



def whichGroup(_pure_alname, _df_group):

    """
    used in "MakeGroupedAllelelist()" as a submodule.

    Allele_group table이 dictionary형태로 _df_group으로 주어지고, 어떤 Allele이 주어졌을때(ex. "01:01:01:01", "A*"는 떼고),
    얘가 어떤 group에 속하는지를 알려줌.

    """

    if not isinstance(_df_group, pd.DataFrame):
        print("[whichGroup]: given _df_group is not DataFrame class.")
        sys.exit()


    flag_group = _df_group.loc[:, "alleles"].apply(lambda x : True if _pure_alname in x else False)

    if flag_group.any():

        # 찾아온게 적어도 하나라도 있으면
        Found_group = _df_group.loc[flag_group]

        if len(Found_group) == 1:
            # 정상적이라면 찾아온 group이 한개여야 함.

            group = Found_group.iat[0, 1] if bool(Found_group.iat[0, 1]) else Found_group.iat[0, 0].pop()
            return group

        elif len(Found_group) > 1:
            # Thankfully, there was no case like this.
            return -2

    else:
        # return -1
        return 0
        # (2018. 7. 1) modified to return '0' instead of '-1' if found nothing.




if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################
    
        ClassifyGroups.py
        
        This script converts a file with allels group information to utilizable file.
        
        It takes (1) "Allelelist.xxxx.txt", (2) "hla_nom_g.txt", (3) "hla_nom_p.txt", (4) IMGT-HLA version,
        (5) "HLA_DICTIONARY_AA", (6) "HLA_DICTIONARY_SNPS" (total 6 things) as an inputs, and generated 
        "Integrated Allele Table" to be used in "NomenCleaner.py" module.
        
        
        made by Wanson Choi.
    
    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)

    parser.add_argument("-al", help="\n\"Allelelist.txt\" file.\n\n", required=True)
    parser.add_argument("-Ggroup", help="\n\"hla_nom_g.txt\" file.\n\n", required=True)
    parser.add_argument("-Pgroup", help="\n\"hla_nom_p.txt\" file.\n\n", required=True)
    parser.add_argument("-imgt", help="\nThe version that you are working with \"Allelelist.txt\", \"hla_nom_g.txt\", \"hla_nom_p.txt\", etc.\n\n",
                        required=True)

    parser.add_argument("-dict-AA", help="\nHLA Sequene information for Amino Acids(ex. \"HLA_DICTIONARY_AA.hg18.imgt370.txt\").\n\n", required=True)
    parser.add_argument("-dict-SNPS", help="\nHLA Sequence information for SNPS(ex. \"HLA_DICTIONARY_SNPS.hg18.imgt370.txt\").\n\n", required=True)




    # <for Publish>
    args = parser.parse_args()

    # <for Testing>

    # hla_nom_p.txt (imgt3320)
    # args = parser.parse_args(["-al", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/Allelelist.3320.txt",
    #                           "-Ggroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/hla_nom_g.txt",
    #                           "-Pgroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/hla_nom_p.txt",
    #                           "-o", "../INTEGRATED_ALLELE_TABLE",
    #                           "-imgt", "3320"])

    # (2018. 7. 12.) imgt370 test
    # args = parser.parse_args(["-al", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/IMGTHLA370/Allelelist.370.txt",
    #                           "-Ggroup", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/IMGTHLA370/wmda/hla_nom_g.txt",
    #                           "-Pgroup", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/IMGTHLA370/wmda/hla_nom_p.txt",
    #                           "-o", "./INTEGRATED_ALLELE_TABLE",
    #                           "-imgt", "370"])

    # # (2018. 7. 15.) "HLA_DICTIONARY" files are introduced to be used as inputs.
    # args = parser.parse_args(["-al", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/IMGTHLA370/Allelelist.370.txt",
    #                           "-Ggroup", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/IMGTHLA370/wmda/hla_nom_g.txt",
    #                           "-Pgroup", "/Users/wansun/Git_Projects/MakeDictionary_RECODE/makedictionary/data/IMGTHLA370/wmda/hla_nom_p.txt",
    #                           "-dict-AA", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-dict-SNPS", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-o", "../INTEGRATED_ALLELE_TABLE",
    #                           "-imgt", "370"])



    print(args)


    ClassifyGroups(args.al, args.Ggroup, args.Pgroup, args.dict_AA, args.dict_SNPS, args.o, args.imgt)




# ========== < deprecated > ==========



# def WhatYangMade(_p_input, _group, _out):
#
#     """
#
#     This function generates an allele information table which Yang has in her project("G-group_alleles.txt").
#
#     Output will be either
#
#     [1] only groups which have at least 2 elements per goup.
#     A*01:01:01G 01:01:01:01/01:01:01:02N/01:01:01:03/01:01:01:04/01:01:01:05/01:01:01:06/01:01:01:07/01:01:01:08/01:01:01:09/01:01:01:10/01:01:38L/01:01:51/01:04N/01:22N/01:32/01:37/01:45/01:56N/01:81/01:87N/01:103/01:107/01:109/01:132/01:141/01:142/01:155/01:177/01:212/01:217/01:234
#     A*01:03:01G 01:03:01:01/01:03:01:02
#     A*02:01:01G 02:01:01:01/02:01:01:02L/02:01:01:03/02:01:01:04/02:01:01:05/02:01:01:06/02:01:01:07/02:01:01:08/02:01:01:09/02:01:01:10/02:01:01:11/02:01:01:12/02:01:01:13/02:01:01:14/02:01:01:15/02:01:01:16/02:01:01:17/02:01:01:18/02:01:01:19/02:01:01:20/02:01:01:21/02:01:01:22/02:01:01:23/02:01:01:24/02:01:01:25/02:01:01:26/02:01:08/02:01:11/02:01:14Q/02:01:15/02:01:21/02:01:48/02:01:50/02:01:79/02:01:80/02:01:89/02:01:97/02:01:98/02:01:99/02:01:104/02:01:130/02:09/02:43N/02:66/02:75/02:83N/02:89/02:97:01/02:97:02/02:132/02:134/02:140/02:241/02:252/02:256/02:266/02:291/02:294/02:305N/02:327/02:329/02:356N/02:357/02:397/02:411/02:446/02:455/02:469/02:481/02:538/02:559/02:607/02:608N/02:614/02:629/02:642/02:665/02:675N
#     A*02:02:01G 02:02:01:01/02:02:01:02
#     A*02:03:01G 02:03:01/02:03:06/02:253/02:264/02:370/02:480/02:505/02:557
#     A*02:04:01G 02:04/02:664
#     ...
#
#     or
#
#     [2] whole groups no matter how many elements in each group.
#     A*01:01:01G	01:01:01:01/01:01:01:02N/01:01:01:03/01:01:01:04/01:01:01:05/01:01:01:06/01:01:01:07/01:01:01:08/01:01:01:09/01:01:01:10/01:01:01:11/01:01:01:12/01:01:01:13/01:01:38L/01:01:51/01:01:83/01:01:84/01:04N/01:22N/01:32/01:37/01:45/01:56N/01:81/01:87N/01:103/01:107/01:109/01:132/01:141/01:142/01:155/01:177/01:212/01:217/01:234/01:237/01:246/01:248Q/01:249/01:251/01:252/01:253
#     A*01:01:02	01:01:02
#     A*01:01:03	01:01:03
#     A*01:01:04	01:01:04
#     A*01:01:05	01:01:05
#     A*01:01:06	01:01:06
#     ...
#
#     Also, "G(3-fields)" or "P(2-fields)" will be determined by a file given to argument `_group`.
#
#     There will be no arguemnt for imgt version(ex. 3320) because it will be determined by header information of
#     given "Allelelist.***.txt" file.
#
#     """
#
#
#     ########## < Core Variables > ##########
#
#     HLA_names  = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
#
#     # General Information
#     ONLY_GROUPS = None
#     WHOLE_GROUPS = None
#
#     # for Utilization
#     ONLY_GROUPS_formal = None
#     WHOLE_GROUPS_formal = None
#
#     std_MAIN_PROCESS_NAME = "[%s]: " % (os.path.basename(__file__))
#
#     print(std_MAIN_PROCESS_NAME+"Reading \"hla_nom_{0}.txt\" file.\n\n".format(_group))
#
#
#     ########## < Reading "hla_nom_g.txt" or "hla_nom_p.txt" file. > ##########
#
#     f = open(_p_input, 'r')
#     lines = f.readlines()
#
#     ### Header info.
#
#     # lines[0:6]  first 6 lines are header information.
#     # Among them, we will take only `data` and `version` info.
#
#     date = lines[1]
#     version = lines[2].rstrip('\n').split(' ')[-1]
#     version = re.sub(pattern='\.', string=version, repl='')
#
#     # print("[ClassifyGroups]: Header information ")
#     print("\nDate : {0}\nVersion : {1}\n".format(date, version))
#
#     sr_hla_g = pd.Series(lines[6:]).apply(lambda x: x.rstrip('\n')).apply(lambda x: x.split(';'))
#
#     df_hla_g = pd.DataFrame(sr_hla_g.tolist(), columns=["HLA", "Alleles", "Group"])
#
#     print("\nPrimitive DataFrame for Grouped allele information.\n")
#     print(df_hla_g.head())
#
#     ########## < Making `ONLY_GROUPS` DataFrame. > ##########
#
#     f_vector = df_hla_g.loc[:, "Group"].apply(lambda x: True if bool(x) else False)
#
#     ONLY_GROUPS = df_hla_g.loc[f_vector, :]
#
#     print("\nDataFrame for `ONLY_GROUPS`.\n")
#     print(ONLY_GROUPS.head())
#
#     col_temp = ONLY_GROUPS.loc[:, ["HLA", "Group"]].apply(''.join, axis=1)
#
#     print("\ncol_temp\n")
#     print(col_temp.head())
#
#     ONLY_GROUPS = pd.concat([col_temp, ONLY_GROUPS.loc[:, "Alleles"]], axis=1)
#     print("\n`ONLY_GROUPS` done.\n")
#     print(ONLY_GROUPS.head())
#
#     ### Exporting `ONLY_GROUPS`
#     ONLY_GROUPS.to_csv('.'.join([_out, "imgt" + version, "alleles.{0}_group.txt".format(_group)]), sep='\t',
#                        header=False, index=False)
#
#     ########## < Making `WHOLE_GROUPS` DataFrame. > ##########
#
#     col_temp = []
#
#     for i in range(0, len(f_vector)):
#
#         if f_vector.iat[i]:
#             # Grouped인 경우
#             col_temp.append(df_hla_g.iat[i, 2])
#         else:
#             # Singleton 인 경우
#             col_temp.append(df_hla_g.iat[i, 1])
#
#     print("\nChoosing group name.\n")
#     print(col_temp)
#     col_temp = pd.concat([df_hla_g.loc[:, "HLA"], pd.Series(col_temp)], axis=1).apply(''.join, axis=1)
#     print(col_temp.head())
#
#     WHOLE_GROUPS = pd.concat([col_temp, df_hla_g.loc[:, "Alleles"]], axis=1)
#     print("\n`WHOLE_GROUPS` done.\n")
#     print(WHOLE_GROUPS.head())
#
#     ### Exporting `WHOLE_GROUPS`
#     WHOLE_GROUPS.to_csv('.'.join([_out, "imgt" + version, "whole.{0}_group.txt".format(_group)]), sep='\t',
#                         header=False, index=False)
#
#     return 0




# def MakeGroupedAllelelist(_p_allelelist, _p_Ggroup, _p_Pgroup, _out):
#
#     """
#
#     :param _p_allelelist: "Allelelist.3320.txt".
#     :param _p_Ggroup: "hla_nom_g.txt".
#     :param _p_Pgroup: "hla_nom_p.txt".
#     :param _out: output file name prefix.
#     :return:
#     """
#
#
#
#     ########## < Core Variables > ##########
#
#     HLA_names  = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
#
#
#
#     ########## < Reading "Allelelist.txt" file > ##########
#
#     f = open(_p_allelelist, 'r')
#     lines = f.readlines()
#
#     header = lines[6].rstrip('\n').split(',')
#     major_lines = lines[7:]
#
#
#
#     df_Allelelist = pd.DataFrame([item.rstrip('\n').split(',') for item in major_lines], columns=header)
#
#     print("\nLoaded Allelelist file.\n")
#     print(df_Allelelist.head())
#
#
#
#     ########## < Reading "hla_nom_g.txt" file > ##########
#
#     f = open(_p_Ggroup, 'r')
#     lines = f.readlines()
#     major_lines = lines[6:]
#
#     sr_hla_g = pd.Series(major_lines).apply(lambda x: x.rstrip('\n')).apply(lambda x: x.split(';'))
#     df_hla_g = pd.DataFrame(sr_hla_g.tolist(), columns=["HLA", "alleles", "Group"])
#
#     # column trimming
#     df_hla_g.loc[:, "HLA"] = df_hla_g.loc[:, "HLA"].apply(lambda x : x.rstrip('*'))
#     df_hla_g.loc[:, "alleles"] = df_hla_g.loc[:, "alleles"].apply(lambda x : x.split('/'))
#
#     # setting index with "HLA" column values
#     df_hla_g = df_hla_g.set_index("HLA")
#
#     print("\nLoaded g group file.\n")
#     print(df_hla_g.head())
#
#
#
#     dict_hla_g = {HLA_names[i] : df_hla_g.loc[HLA_names[i], :] for i in range(0, len(HLA_names))}
#
#
#
#     ########## < Reading "hla_nom_p.txt" file > ##########
#
#     f = open(_p_Pgroup, 'r')
#     lines = f.readlines()
#     major_lines = lines[6:]
#
#     sr_hla_p = pd.Series(major_lines).apply(lambda x: x.rstrip('\n')).apply(lambda x: x.split(';'))
#     df_hla_p = pd.DataFrame(sr_hla_p.tolist(), columns=["HLA", "alleles", "Group"])
#
#     # column trimming
#     df_hla_p.loc[:, "HLA"] = df_hla_p.loc[:, "HLA"].apply(lambda x : x.rstrip('*'))
#     df_hla_p.loc[:, "alleles"] = df_hla_p.loc[:, "alleles"].apply(lambda x : x.split('/'))
#
#     # setting index with "HLA" column values
#     df_hla_p = df_hla_p.set_index("HLA")
#
#
#     print("\nLoaded p group file.\n")
#     print(df_hla_p.head(10))
#
#
#
#     dict_hla_p = {HLA_names[i] : df_hla_p.loc[HLA_names[i], :] for i in range(0, len(HLA_names))}
#
#
#
#
#     ########## < Matching Job > ##########
#
#     # 이제 df_Allelelist의 "alleles" 컬럼의 원소에 대해 for문을 돌거임.
#     # 그런데 그전에 우선 HLA_names 애들만 건저내야하 할듯(ex. "W*01:01" 이런애들)
#
#     splitted_alleles = pd.DataFrame(df_Allelelist.loc[:, "Allele"].apply(lambda x : x.split('*')).tolist(), columns=["HLA", "Pure_Alleles"])
#     # print(splitted_alleles.head())
#
#     df_Allelelist_filtered = pd.concat([df_Allelelist.loc[:, "Allele"], splitted_alleles], axis=1).set_index("HLA").loc[HLA_names, :] # Filtering only major HLA genes.
#     df_Allelelist_filtered = df_Allelelist_filtered.reset_index("HLA").set_index("Allele")
#     print(df_Allelelist_filtered.head())
#
#
#     Found_Gs = pd.Series([whichGroup(df_Allelelist_filtered.iat[i, 1], dict_hla_g[df_Allelelist_filtered.iat[i, 0]]) for i in range(0, len(df_Allelelist_filtered))])
#     print(Found_Gs)
#     # Found_Gs.to_csv("./Found_Gs.txt")
#
#     Found_Ps = pd.Series([whichGroup(df_Allelelist_filtered.iat[i, 1], dict_hla_p[df_Allelelist_filtered.iat[i, 0]]) for i in range(0, len(df_Allelelist_filtered))])
#     print(Found_Ps)
#     # Found_Ps.to_csv("./Found_Ps.txt")
#
#
#
#     df_OUTPUT = pd.concat([Found_Gs, Found_Ps], axis=1)
#     df_OUTPUT.index = df_Allelelist_filtered.index
#     df_OUTPUT.columns = ["G_group", "P_group"]
#     df_OUTPUT.to_csv('./Allele_Group_table_1st.txt', sep='\t', header=True, index=True)
#
#     print("\nOutput\n")
#     print(df_OUTPUT.head(50))
#
#
#
#
#     # # for i in range(0, len(alleles))
#     # for i in range(0, 30):
#     #
#     #     hla, pure_alname = alleles[i].split('*')
#     #
#     #     print(whichGroup(pure_alname, dict_hla_g[hla]))
#     #     # whichGroup(pure_alname, dict_hla_p[hla])
#
#
#     return 0