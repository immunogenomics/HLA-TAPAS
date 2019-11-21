# -*- coding: utf-8 -*-


import os, sys, re
import argparse, textwrap
import pandas as pd
# from collections import OrderedDict

def NomenCleaner(_p_ped, _ped_descriptor, _p_iat, _out, _field_format, _f_NoCaption = False):
    """
    NomenCleaner.py

    Taking (1) "*.ped" file with HLA alleles and (2) Integrated_Allele_Table("*.iat") file,
    this module generates ped file with given output format which user selected(1,2,3,4-field or P,G group).

    """

    ########## < Core Variables > ##########

    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    header_ped = ["FamID", "IdivID", "P_ID", "M_ID", "Sex", "Phe"]

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME + "Init.")

    PED = None
    IAT = None

    ########## < Loading "*.ped" file > ##########

    PED = pd.read_table(_p_ped, sep='\t', header=None, dtype=str,
                        names=header_ped + [item + "_" + str(i) for item in HLA_names for i in range(1, 3)]).set_index((header_ped))

    print(std_MAIN_PROCESS_NAME + "Loaded \"*.ped\" file.")
    print(PED.head())

    ########## < Loading "*.iat" file > ##########

    IAT = pd.read_table(_p_iat, sep='\t', header=0, dtype=str)

    print(std_MAIN_PROCESS_NAME + "Loaded \"*.iat\" file.\n")
    print(IAT.head())

    IAT = pd.concat([pd.DataFrame(IAT.loc[:, "Allele"].apply(lambda x: x.split('*')).tolist(), columns=["HLA", "Allele"]), IAT.loc[:, ["G_group", "P_group"]]], axis=1).set_index("HLA")
    print("\nNew IAT\n")
    print(IAT.head())

    IAT_dict = {HLA_names[i]: IAT.loc[HLA_names[i], :].set_index("Allele") for i in range(0, len(HLA_names))}

    print("\nIAT divided by HLA gene names.\n")
    print(IAT_dict["C"].head(10))

    """
    (2018. 7. 2.)
    As i decided to classify taking input "*.ped" file into 3 categories, (1) Standard 4-field allele ped file,
    (2) G-group allele ped file, (3) P-group allele ped file, the main process in this module also will be divided 
    to apply those changes.

    The variable `_p_ped_descriptor` represents which type of ped file is given among those 3 options. It has a value
    1,2 or 3 which consecutively represents type (1), (2), or (3) mentioned above. 

    I will divide the main process into two major process which basedon type (1) and type (2), (3).
    When P or G group allele ped file is given, we can assum that those allele didn't go through trimming job.
    At least I will define this module as a function that takes P or G group allele which didn't go through trimming job.
    So, In case P or G group allele is given, by just checking "Gene caption or not" and "has Double-colon or not",
    we can find any format of allele with `Integrated_Allele_Table(*.iat)` file. I will divide main process into two,
    and first part will have the process that i originally made and second part will have the new code that deals with
    what i mentioned now(The new code that works taking advantage of the fact that P or G group allele didn't go through trimming job).   
    """

    ########## < Checking and Transforming each alleles(Main For loop) > ##########

    l_temp = []

    for i in range(0, len(HLA_names)):  # iterating over `HLA_names`
        # for i in range(0, 1): # iterating over `HLA_names`

        curr_hla_name = HLA_names[i]
        curr_IAT_dict = IAT_dict[curr_hla_name]

        # # set of alleles in "Allelelist.txt"
        # idx_list = curr_IAT_dict.index.tolist()

        l_temp.append(PED.loc[:, [curr_hla_name + "_1", curr_hla_name + "_2"]].applymap(
            lambda x: Main_Transformation(x, curr_hla_name, curr_IAT_dict, _ped_descriptor) if x != "0" else "0"))

    df_TRANSFORMED_4field = pd.concat(l_temp, axis=1)

    print(std_MAIN_PROCESS_NAME + "Transformed ped DataFrame.\n")
    print(df_TRANSFORMED_4field.head())

    ### [4] Choosing output format based on `_field_format`

    if _field_format == 1:

        print(std_MAIN_PROCESS_NAME + "Output as 1-field.\n")

        p = re.compile("\d{2,3}[A-Z]?")
        df_TRANSFORMED_1field = df_TRANSFORMED_4field.applymap(
            lambda x: p.match(string=x).group() if bool(p.match(string=x)) else x)

        if not _f_NoCaption:
            df_TRANSFORMED_1field = pd.concat([df_TRANSFORMED_1field.iloc[:, [2*i, 2*i+1]].applymap(lambda x : '*'.join([HLA_names[i], x]) if x != "0" else x) for i in range(0, len(HLA_names))], axis=1)
            df_TRANSFORMED_1field.index = df_TRANSFORMED_4field.index

        print(std_MAIN_PROCESS_NAME + "1-field output ped file.\n")
        print(df_TRANSFORMED_1field.head())

        df_TRANSFORMED_1field.to_csv(_out + ".1field.ped", sep='\t', header=False, index=True)


    elif _field_format == 2:

        print(std_MAIN_PROCESS_NAME + "Output as 2-field.\n")

        p = re.compile("\d{2,3}\:\d{2,3}[A-Z]?")
        df_TRANSFORMED_2field = df_TRANSFORMED_4field.applymap(
            lambda x: p.match(string=x).group() if bool(p.match(string=x)) else x)

        if not _f_NoCaption:
            df_TRANSFORMED_2field = pd.concat([df_TRANSFORMED_2field.iloc[:, [2*i, 2*i+1]].applymap(lambda x : '*'.join([HLA_names[i], x]) if x != "0" else x) for i in range(0, len(HLA_names))], axis=1)
            df_TRANSFORMED_2field.index = df_TRANSFORMED_4field.index

        print(std_MAIN_PROCESS_NAME + "2-field output ped file.\n")
        print(df_TRANSFORMED_2field.head())

        df_TRANSFORMED_2field.to_csv(_out + ".2field.ped", sep='\t', header=False, index=True)


    elif _field_format == 3:

        print(std_MAIN_PROCESS_NAME + "Output as 3-field.\n")

        p = re.compile("\d{2,3}\:\d{2,3}\:\d{2,3}[A-Z]?")
        df_TRANSFORMED_3field = df_TRANSFORMED_4field.applymap(
            lambda x: p.match(string=x).group() if bool(p.match(string=x)) else x)

        if not _f_NoCaption:
            df_TRANSFORMED_3field = pd.concat([df_TRANSFORMED_3field.iloc[:, [2*i, 2*i+1]].applymap(lambda x : '*'.join([HLA_names[i], x]) if x != "0" else x) for i in range(0, len(HLA_names))], axis=1)
            df_TRANSFORMED_3field.index = df_TRANSFORMED_4field.index

        print(std_MAIN_PROCESS_NAME + "3-field output ped file.\n")
        print(df_TRANSFORMED_3field.head())

        df_TRANSFORMED_3field.to_csv(_out + ".3field.ped", sep='\t', header=False, index=True)


    elif _field_format == 4:

        print(std_MAIN_PROCESS_NAME + "Output as standard 4-field.\n")

        if not _f_NoCaption:

            df_idx = df_TRANSFORMED_4field.index

            df_TRANSFORMED_4field = pd.concat([df_TRANSFORMED_4field.iloc[:, [2*i, 2*i+1]].applymap(lambda x : '*'.join([HLA_names[i], x]) if x != "0" else x) for i in range(0, len(HLA_names))], axis=1)
            df_TRANSFORMED_4field.index = df_idx

        print(std_MAIN_PROCESS_NAME + "4-field output ped file.\n")
        print(df_TRANSFORMED_4field.head())

        df_TRANSFORMED_4field.to_csv(_out + ".4field.ped", sep='\t', header=False, index=True)


    elif _field_format == 5:

        print(std_MAIN_PROCESS_NAME + "Output as G-group.\n")

        # for i in range(0, len(HLA_names)):
        #     print(df_TRANSFORMED_4field.loc[:, [HLA_names[i] + "_1", HLA_names[i] + "_2"]].)

        df_TRANSFORMED_Ggroup = pd.concat([df_TRANSFORMED_4field.loc[:,
                                           [HLA_names[i] + "_1", HLA_names[i] + "_2"]].applymap(
            lambda x: IAT_dict[HLA_names[i]].loc[x, "G_group"] if str(x) != "0" else x) for i in
                                           range(0, len(HLA_names))], axis=1)

        if not _f_NoCaption:
            df_TRANSFORMED_Ggroup = pd.concat([df_TRANSFORMED_Ggroup.iloc[:, [2*i, 2*i+1]].applymap(lambda x : '*'.join([HLA_names[i], x]) if x != "0" else x) for i in range(0, len(HLA_names))], axis=1)
            df_TRANSFORMED_Ggroup.index = df_TRANSFORMED_4field.index

        print(std_MAIN_PROCESS_NAME + "G-group output ped file.\n")
        print(df_TRANSFORMED_Ggroup.head())

        df_TRANSFORMED_Ggroup.to_csv(_out + ".Ggroup.ped", sep='\t', header=False, index=True)


    elif _field_format == 6:

        print(std_MAIN_PROCESS_NAME + "Output as P-group.\n")

        df_TRANSFORMED_Pgroup = pd.concat([df_TRANSFORMED_4field.loc[:,
                                           [HLA_names[i] + "_1", HLA_names[i] + "_2"]].applymap(
            lambda x: IAT_dict[HLA_names[i]].loc[x, "P_group"] if str(x) != "0" else x) for i in
                                           range(0, len(HLA_names))], axis=1)

        if not _f_NoCaption:
            df_TRANSFORMED_Pgroup = pd.concat([df_TRANSFORMED_Pgroup.iloc[:, [2*i, 2*i+1]].applymap(lambda x : '*'.join([HLA_names[i], x]) if x != "0" else x) for i in range(0, len(HLA_names))], axis=1)
            df_TRANSFORMED_Pgroup.index = df_TRANSFORMED_4field.index

        print(std_MAIN_PROCESS_NAME + "P-group output ped file.\n")
        print(df_TRANSFORMED_Pgroup.head())

        df_TRANSFORMED_Pgroup.to_csv(_out + ".Pgroup.ped", sep='\t', header=False, index=True)

    return 0


def isCaptioned(_the_allele, _hla_name):
    """
    Check whether an allele given as argument `_the_allele` has gene caption(ex. "A*") or not.
    If it does, then return the allele of which gene caption is trimmed off.
    """

    if bool(re.search(string=_the_allele, pattern=_hla_name + "\*")):
        return re.sub(pattern=_hla_name + "\*", string=_the_allele, repl='')
    else:
        return _the_allele


def hasDoubleColon(_the_allele):
    return True if _the_allele.find(':') >= 0 else False


def CHECK_DIGITS(_hla_name, _the_allele, _IAT_Allelelist):
    """

    Perform this function assuming given `_the_allele` neither is captioned or has double-colon.

    (2018. 6. 29.)
    More classification conditions are added.

    (2018. 7. 23.)
    I found the way to deal with a Tag has not been introduced yet.
    Now it is introduced.

    """

    # check whether some single character is tagged along with the given allele.
    p_tag = re.compile(r'.+[A-Z]$')
    hasTag = p_tag.match(_the_allele)

    if hasTag:
        t_name = _the_allele[:-1]
        t_name_tag = _the_allele[-1]
    else:
        t_name = _the_allele
        t_name_tag = -1

    sr_IAT_Allelelist = pd.Series(_IAT_Allelelist)


    if len(t_name) == 2:
        # 1-field / 2-digits
        return t_name[0:2] + (t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 3:
        # 1-field / 3-digits
        return t_name[0:3] + (t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 4:
        # 2 + 2 (+C) = 4

        return ':'.join([t_name[0:2], t_name[2:4]]) + (t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 5:
        # (1) 2 + 3 (+C) = 5
        # (2) 3 + 2 = 5

        p = re.compile('\:'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any = sr_IAT_Allelelist.str.match(p).any()

        return (':'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
            (':'.join([t_name[0:3], t_name[3:5]]) + (t_name_tag if t_name_tag != -1 else ""))

    elif len(t_name) == 6:
        # (1) 2 + 2 + 2 (+C)
        # (2) 3 + 3 (+C)

        p = re.compile('\:'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any = sr_IAT_Allelelist.str.match(p).any()

        return (':'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
            (':'.join([t_name[0:3], t_name[3:6]]) + (t_name_tag if t_name_tag != -1 else ""))

    elif len(t_name) == 7:
        # (1) 2 + 3 + 2 (+C) = 7
        # (2) 2 + 2 + 3 (+C) = 7
        # (3) 3 + 2 + 2 (+C) = 7

        p1 = re.compile('\:'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any1 = sr_IAT_Allelelist.str.match(p1).any()

        if not Found_Any1:
            # Not the case of "(1) 2 + 3 + 2 = 7"
            p2 = re.compile('\:'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any2 = sr_IAT_Allelelist.str.match(p2).any()

            if not Found_Any2:
                # Not the case of "(2) 2 + 2 + 3 = 7"
                p3 = re.compile('\:'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
                Found_Any3 = sr_IAT_Allelelist.str.match(p3).any()

                if not Found_Any3:
                    # Not the case of "(3) 3 + 2 + 2 = 7"
                    # Not Found
                    return "-1"
                else:
                    # The case of "(3) 3 + 2 + 2 = 7"
                    return ':'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else "")
            else:
                # The case of "(2) 2 + 2 + 3 = 7"
                return ':'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else "")
        else:
            # The case of "(1) 2 + 3 + 2 = 7"
            return ':'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else "")


        # return ':'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) if Found_Any else \
        #     ':'.join([t_name[0:3], t_name[3:5], t_name[5:7]])

    elif len(t_name) == 8:
        # (1) 2 + 2 + 2 + 2
        # (2) 3 + 2 + 3
        # (3) 3 + 3 + 2

        p1 = re.compile('\:'.join([t_name[0:2], t_name[2:4], t_name[4:6], t_name[6:8]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any1 = sr_IAT_Allelelist.str.match(p1).any()

        if not Found_Any1:
            # Not the case of "(1) 2 + 2 + 2 + 2"
            p2 = re.compile('\:'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any2 = sr_IAT_Allelelist.str.match(p2).any()

            if not Found_Any2:
                # Not the case of "(2) 3 + 2 + 3"
                p3 = re.compile('\:'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else ""))
                Found_Any3 = sr_IAT_Allelelist.str.match(p3).any()

                if not Found_Any3:
                    return str(-1)

                else:
                    return ':'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else "")
            else:
                return ':'.join([t_name[0:3], t_name[3:5], t_name[5:8]]) + (t_name_tag if t_name_tag != -1 else "")
        else:
            return ':'.join([t_name[0:2], t_name[2:4], t_name[4:6], t_name[6:8]]) + (t_name_tag if t_name_tag != -1 else "")

    elif len(t_name) == 9:
        # (1) 2 + 3 + 2 + 2
        # (2) 3 + 2 + 2 + 2

        p = re.compile('\:'.join([t_name[0:2], t_name[2:5], t_name[5:7], t_name[7:9]]) + (t_name_tag if t_name_tag != -1 else ""))
        Found_Any = sr_IAT_Allelelist.str.match(p).any()

        return (':'.join([t_name[0:2], t_name[2:5], t_name[5:7], t_name[7:9]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
            (':'.join([t_name[0:3], t_name[3:5], t_name[5:7], t_name[7:9]]) + (t_name_tag if t_name_tag != -1 else ""))

    elif len(t_name) == 10:
        # (1) 3 + 3 + 2 + 2

        return ':'.join([t_name[0:3], t_name[3:6], t_name[6:8], t_name[8:10]]) + (t_name_tag if t_name_tag != -1 else "")


def CHECK_DIGITS_PorGgroup(_hla_name, _the_allele, _IAT_Allelelist, _ped_descriptor):
    """

    ### < G-group > ###

    # 2 + 2 (+C) = 4
    - 46:51Q

    # 2 + 3 (+C) = 5
    - 08:124, 07:534
    - 01:247N, 01:250N

    # 3 + 2 = 5
    - 255:01, 527:01, 632:01, 137:01, 338:01

    # 2 + 2 + 2 (+C) = 6
    - 01:01:02
    - 01:01:01G

    # 2 + 3 + 2 = 7
    - 01:146:01
    - 06:127:01G

    # 2 + 2 + 3 = 7
    - 02:01:100, 02:01:101

    # 3 + 2 + 2 = 7
    - 155:01:01


    ### < P-group > ###

    # 2 + 2 = 4
    - 10:24, 11:04
    - 01:01P, 02:01P

    # 2 + 3 = 5
    - 14:110, 14:139
    - 01:146P, 02:101P, 02:610P

    # 3 + 2 = 5
    - 100:01, 119:01
    - 155:01P, 279:01P

    # 2 + 2 + 2 = 6
    - 30:12:01 (2018. 7. 5) I found this case lol!


    """

    # _IAT_Allelelist = pd.Series(_IAT_Allelelist)

    # check whether some single character is tagged along with the given allele.
    p_tag = re.compile(r'.+[A-Z]$')
    hasTag = p_tag.match(_the_allele)

    if hasTag:
        t_name = _the_allele[:-1]
        t_name_tag = _the_allele[-1]
    else:
        t_name = _the_allele
        t_name_tag = -1

    # print("\n[CHECK_DIGITS_PorG] : {0} and {1}".format(t_name, t_name_tag))

    ### In case of "G-group"
    if _ped_descriptor == 2:

        # Main matching job

        if len(t_name) == 4:
            # (1) 2 + 2 (+C) = 4

            return ':'.join([t_name[0:2], t_name[2:4]]) + (t_name_tag if t_name_tag != -1 else "")

        elif len(t_name) == 5:
            # (1) 2 + 3 (+C) = 5
            # (2) 3 + 2 = 5

            p = re.compile('\:'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any = _IAT_Allelelist.str.match(p).any()

            return (':'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else \
                (':'.join([t_name[0:3], t_name[3:5]]) + (t_name_tag if t_name_tag != -1 else ""))

        elif len(t_name) == 6:
            # (1) 2 + 2 + 2 (+C) = 6

            return ':'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (t_name_tag if t_name_tag != -1 else "")

        elif len(t_name) == 7:
            # (1) 2 + 3 + 2 = 7
            # (2) 2 + 2 + 3 = 7
            # (3) 3 + 2 + 2 = 7

            p1 = re.compile('\:'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any1 = _IAT_Allelelist.str.match(p1).any()

            if not Found_Any1:
                # Not the case of "(1) 2 + 3 + 2 = 7"
                p2 = re.compile('\:'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else ""))
                Found_Any2 = _IAT_Allelelist.str.match(p2).any()

                if not Found_Any2:
                    # Not the case of "(2) 2 + 2 + 3 = 7"
                    p3 = re.compile('\:'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else ""))
                    Found_Any3 = _IAT_Allelelist.str.match(p3).any()

                    if not Found_Any3:
                        return "-1"

                    else:
                        return ':'.join([t_name[0:3], t_name[3:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else "")
                else:
                    return ':'.join([t_name[0:2], t_name[2:4], t_name[4:7]]) + (t_name_tag if t_name_tag != -1 else "")
            else:
                return ':'.join([t_name[0:2], t_name[2:5], t_name[5:7]]) + (t_name_tag if t_name_tag != -1 else "")

        else:
            print("\nNo match for G_group.\n")
            return "-1"


    ### In case of "P-group"
    elif _ped_descriptor == 3:

        # Main Matching Job

        if len(t_name) == 4:
            # 2 + 2 = 4
            return ':'.join([t_name[0:2], t_name[2:4]]) + (t_name_tag if t_name_tag != -1 else "")

        elif len(t_name) == 5:
            # 2 + 3 = 5
            # 3 + 2 = 5

            p = re.compile('\:'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else ""))
            Found_Any = _IAT_Allelelist.str.match(p).any()

            return (':'.join([t_name[0:2], t_name[2:5]]) + (t_name_tag if t_name_tag != -1 else "")) if Found_Any else (
                    ':'.join([t_name[0:3], t_name[3:5]]) + (t_name_tag if t_name_tag != -1 else ""))

        elif len(t_name) == 6:
            # 2 + 2 + 2
            # 3 + 3 ... Not Found.

            return ':'.join([t_name[0:2], t_name[2:4], t_name[4:6]]) + (t_name_tag if t_name_tag != -1 else "")


        else:
            print("\nNo match for P_group.\n")
            return "-1"


def Find_1st_Allele(_the_allele, _IAT_Allelelist):
    """

    In case where given allele has double-colon(implying its digits are determined) but not in a form of complete allele name in "Allelelist.txt" file,
    then i will classify this allele as "Non-deterministic 4-field allele".

    So, that incomplete allele("Non-deterministic 4-field allele") will be matched by "re.match" function and return 1st element.
    """

    if not (isinstance(_IAT_Allelelist, list) or isinstance(_IAT_Allelelist, tuple)):
        print("Given _IAT_Allelelist is not a list or tuple.\n")
        sys.exit()

    t_sr = pd.Series(_IAT_Allelelist)


    if 2 <= len(_the_allele) <= 5:

        """
        (2018. 7. 23.)
        This classification block was introduced as I found some unproper exceptions.
        For example, When the allele "B*1501" is given, the transformed result is "15:170" not "B*15:17:01:01".
        It means digit information isn't fully considered in this function, i.e. this part which makes regular expression pattern.
        
        As a solution, I decided to add more classification related to `len(_the_allele)`.
        """

        Flag_Found = t_sr.str.match(str(_the_allele)+":")

        if Flag_Found.any():
            return t_sr.loc[Flag_Found].iat[0]
        else:
            # One more time
            Flag_Found = t_sr.str.match(str(_the_allele))

            if Flag_Found.any():
                return t_sr.loc[Flag_Found].iat[0]
            else:
                return 0


    else:
        p = re.compile(str(_the_allele))
        # Flag_Found = t_sr.apply(lambda x: bool(p.match(string=x)))
        Flag_Found = t_sr.str.match(str(_the_allele))

        if Flag_Found.any():
            return t_sr.loc[Flag_Found].iat[0]
        else:
            return 0


def Find_1st_Allele_PorGgroup(_the_allele, _df_IAT_Allele, _ped_descriptor):
    """
    If given allele is P or G-group allele, then we can just find corresponding allele using ".loc[]" function.

    """

    if not isinstance(_df_IAT_Allele, pd.DataFrame):
        print("Given `_df_IAT_Allelelist` is not a pandas.DataFrame.\n")
        sys.exit()

    df_temp = _df_IAT_Allele.reset_index().set_index("G_group" if _ped_descriptor == 2 else "P_group")

    Found_df = df_temp.loc[_the_allele, "Allele"]
    # (2018. 7. 4.) Found result `Found_df` could be 'str' object. So, you shouldn't use ".iat" function.

    # return str(Found_df.iat[0]) if len(Found_df) > 0 else "0"

    if isinstance(Found_df, pd.Series):
        # ex. "01:01:01G" => ["01:01:01:01", "01:01:01:02N", "01:01:01:03", ..., "01:253"]
        # Multiple results retured as "Series" object.
        return str(Found_df.iat[0])
    elif isinstance(Found_df, str):
        return Found_df
    else:
        return "0"


def Main_Transformation(_single_allele, _hla_name, _df_IAT_Allelelist, _ped_descriptor):
    """
    Main Processes which were located in main for loop. They are now moved to this function.
    From now on, those main processes will be preformed to single allele.
    """

    _IAT_Allelelist = _df_IAT_Allelelist.index.tolist()

    ### [1] Checking whether given allele is captioned or not.

    _al_filetered1 = isCaptioned(_single_allele, _hla_name)

    ### [2] Checking whether given allele has double-colon or not.

    _al_filetered2 = ""

    if _ped_descriptor == 1:
        # When given ped file is standard 4-field allele(where `_ped_descriptor` == 1).
        _al_filetered2 = CHECK_DIGITS(_hla_name, _al_filetered1, _IAT_Allelelist) if not hasDoubleColon(_al_filetered1) else _al_filetered1

    elif _ped_descriptor == 2:
        # When given ped file is P or G-group allele(where `_ped_descriptor` == 2 or 3).
        _al_filetered2 = CHECK_DIGITS_PorGgroup(_hla_name, _al_filetered1, _df_IAT_Allelelist.loc[:, "G_group"], _ped_descriptor) if not hasDoubleColon(_al_filetered1) else _al_filetered1

    elif _ped_descriptor == 3:
        _al_filetered2 = CHECK_DIGITS_PorGgroup(_hla_name, _al_filetered1, _df_IAT_Allelelist.loc[:, "P_group"], _ped_descriptor) if not hasDoubleColon(_al_filetered1) else _al_filetered1

    ### [3] Digit-Checking to determine proper number of digits per field.

    _al_filetered3 = ""

    if _ped_descriptor < 2:
        _al_filetered3 = Find_1st_Allele(_al_filetered2, _IAT_Allelelist)
    elif _ped_descriptor >= 2:
        _al_filetered3 = Find_1st_Allele_PorGgroup(_al_filetered2, _df_IAT_Allelelist, _ped_descriptor)

    return str(_al_filetered3)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

     NomenCleaner.py


    #########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    # Input (1) : *.ped file
    PED_TYPE = parser.add_mutually_exclusive_group(required=True)
    PED_TYPE.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", dest="ped",
                          default="Not_given")
    PED_TYPE.add_argument("-ped-Ggroup", help="\nHLA Type Data(G-group allele \"*.ped\" file).\n\n", dest="ped_G",
                          default="Not_given")
    PED_TYPE.add_argument("-ped-Pgroup", help="\nHLA Type Data(P-group allele \"*.ped\" file).\n\n", dest="ped_P",
                          default="Not_given")

    # Input (2) : *.iat file
    parser.add_argument("-iat", help="\nIntegrated Allele Table file(*.iat).\n\n", required=True)

    # Ouptut Prefix
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    # Output format
    output_digit_selection = parser.add_mutually_exclusive_group(required=True)
    output_digit_selection.add_argument("--1field", help="\nOutput ped file as '1-field' format.\n\n",
                                        action="store_true", dest="oneF")
    output_digit_selection.add_argument("--2field", help="\nOutput ped file as '2-field' format.\n\n",
                                        action="store_true", dest="twoF")
    output_digit_selection.add_argument("--3field", help="\nOutput ped file as '3-field' format.\n\n",
                                        action="store_true", dest="threeF")
    output_digit_selection.add_argument("--4field",
                                        help="\nOutput ped file as '4-field(Current Standard Names)' format.\n\n",
                                        action="store_true", dest="fourF")
    output_digit_selection.add_argument("--G-group", help="\nOutput ped file as 'G-group' format.\n\n",
                                        action="store_true")
    output_digit_selection.add_argument("--P-group", help="\nOutput ped file as 'P-group' format.\n\n",
                                        action="store_true")

    # Flag to remove HLA gene caption.
    parser.add_argument("--NoCaption", help="\nOutput without HLA gene(ex. \"A*\").\n\n", action='store_true')



    ##### <for Test> #####

    # ==========< Input as standard 4-field allele >==========

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/COATING_PED/HAPMAP_CEU_HLA.5digits.example.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/CCCCOATED_PED_TEST"
    #                           ])

    # # Standard 4-field ped file(500 sample)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/NomenCleaner/DummyPED.standard.500.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/COATING_PED_test_standard.500",
    #                           "--4field"
    #                           ])

    # 2-field ped file(500 sample)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/NomenCleaner/DummyPED.standard.500.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/COATING_PED_test.500.NC",
    #                           "--2field",
    #                           "--NoCaption"
    #                           ])

    # (2018. 7. 9.)
    # args = parser.parse_args(["-ped",
    #                           "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/NomenCleaner/DummyPED.standard.100.ped",
    #                           "-iat",
    #                           "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o",
    #                           "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/Transformed.standard.500",
    #                           "--G-group"
    #                           ])

    # (2018. 7. 12.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt370.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/HAPMAP_CEU_HLA",
    #                           "--4field",
    #                           # "--NoCaption"
    #                           ])

    # ==========< Input as G-group allele >==========

    # # G-group ped file.(No cap, Double-Colon, 100 rows)
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/COATING_PED/DummyPED.Ggroup.Ncap.dc.100.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/CCCCOATED_PED_TEST.Ggroup.Ncap.dc.100",
    #                           "--4field"
    #                           ])

    # # G-group ped file.(No cap, No Double-Colon, 100 rows)
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/COATING_PED/DummyPED.Ggroup.Ncap.Ndc.100.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/CCCCOATED_PED_TEST.Ggroup.Ncap.Ndc.100",
    #                           "--4field"
    #                           ])

    # # G-group ped file.(No cap, Double-Colon, 500 rows)
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/COATING_PED/DummyPED.Ggroup.Ncap.dc.500.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/CCCCOATED_PED_TEST.Ggroup.Ncap.dc.500",
    #                           "--4field"
    #                           ])

    # # G-group ped file.(No cap, No Double-Colon, 500 rows)
    # args = parser.parse_args(["-ped-Ggroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/COATING_PED/DummyPED.Ggroup.Ncap.Ndc.500.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/CCCCOATED_PED_TEST.Ggroup.Ncap.Ndc.500",
    #                           "--4field"
    #                           ])

    # ==========< Input as P-group allele >==========

    # # P-group ped file.(No cap, Double-Colon, 500 rows)
    # args = parser.parse_args(["-ped-Pgroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/COATING_PED/DummyPED.Pgroup.Ncap.dc.500.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/CCCCOATED_PED_TEST.Pgroup.Ncap.dc.500",
    #                           "--4field"
    #                           ])

    # # P-group ped file.(No cap, No Double-Colon, 500 rows)
    # args = parser.parse_args(["-ped-Pgroup", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/COATING_PED/DummyPED.Pgroup.Ncap.Ndc.500.ped",
    #                           "-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/CCCCOATED_PED_TEST.Pgroup.Ncap.Ndc.500",
    #                           "--4field"
    #                           ])

    ##### <for Publication> #####

    args = parser.parse_args()

    ### Additional Argument processing

    ## Output Format Flags
    FILE_FORMAT = 1 if args.oneF else 2 if args.twoF else 3 if args.threeF else 4 if args.fourF else 5 if args.G_group else 6 if args.P_group else -1

    ## Which type of ped file given?
    _p_ped = -1
    _p_ped_descriptor = -1

    if args.ped != "Not_given":
        # Standard 4-field *.ped file given
        _p_ped = args.ped
        _p_ped_descriptor = 1
    elif args.ped_G != "Not_given":
        # G-group *.ped file given
        _p_ped = args.ped_G
        _p_ped_descriptor = 2
    elif args.ped_P != "Not_given":
        # P-group *.ped file given
        _p_ped = args.ped_P
        _p_ped_descriptor = 3
    else:
        # Assuming at least three of them given, there won't be the case which comes to here.
        print("\nArgument for input *.ped file has a problem. Please check it again.\n")
        sys.exit()

    ## If input ped file is given as G-group(or P-group), then user can't choose output format as same G-group(as same to P-group)
    if (_p_ped_descriptor == 2 and FILE_FORMAT == 5):
        print("\nYou just asked pointless transformation(Transformation G-group to G-group is pointless).")
        print("Skip this Transformation Request.\n")
        sys.exit()

        # (2018. 7. 10.) P to P or G to G 여도 NoDoubleColon인 경우는 할 수 있게 해줘야할듯.

    if (_p_ped_descriptor == 3 and FILE_FORMAT == 6):
        print("\nYou just asked pointless transformation(Transformation P-group to P-group is pointless).")
        print("Skip this Transformation Request.\n")
        sys.exit()

    print(args)


    # Implementing Main Function
    NomenCleaner(_p_ped, _p_ped_descriptor, args.iat, args.o, FILE_FORMAT, _f_NoCaption=args.NoCaption)

