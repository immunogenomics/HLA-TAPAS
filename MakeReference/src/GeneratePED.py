# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd
import random

def GeneratePED(_p_iat, _out, _N,
                _f_asStandard = False, _f_asGgroup = False, _f_asPgroup = False,
                _f_caption = -1, _f_trimming = -1, _f_DC = -1):

    """

    """

    ########## < Core Variables > ##########

    HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    header_ped = ["FamID", "IdivID", "P_ID", "M_ID", "Sex", "Phe"]

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME + "Init.\n")

    _N = int(_N) # Type casting `_N` arguemnt from "str" to "int".

    IAT = pd.read_table(_p_iat, sep='\t', header=0, dtype=str)
    IAT = pd.concat([pd.DataFrame(IAT.loc[:, "Allele"].str.split('*').tolist(), columns=["HLA", "Allele"]), IAT.loc[:, ["G_group", "P_group"]]], axis=1).set_index("HLA")

    IAT_dict = {HLA_names[i]: IAT.loc[HLA_names[i], :] for i in range(0, len(HLA_names))}

    print(std_MAIN_PROCESS_NAME + "Loaded \"*.iat\" file.\n")
    print(IAT.head())

    ### Generating Alleles

    l_forChoice = tuple()

    if (_f_asStandard and _f_asGgroup and _f_asPgroup):
        l_forChoice = (0, 1, 2)

    elif _f_asStandard and not _f_asGgroup and not _f_asPgroup:
        # only Standard(4-field) alleles.
        l_forChoice = (0,)
    elif not _f_asStandard and _f_asGgroup and not _f_asPgroup:
        # only G-Group alleles
        l_forChoice = (1,)
    elif not _f_asStandard and not _f_asGgroup and _f_asPgroup:
        # only P-Group alleles
        l_forChoice = (2,)

    elif _f_asStandard and _f_asGgroup and not _f_asPgroup:
        # only Standard(4-field) alleles.
        l_forChoice = (0, 1)
    elif _f_asStandard and not _f_asGgroup and _f_asPgroup:
        # only Standard(4-field) alleles.
        l_forChoice = (0, 2)
    elif not _f_asStandard and _f_asGgroup and _f_asPgroup:
        # only Standard(4-field) alleles.
        l_forChoice = (1, 2)

    elif (not _f_asStandard and not _f_asGgroup and not _f_asPgroup):
        # Most of case will fall into this category
        # only Standard(4-field) alleles.
        l_forChoice = (0,)

    ########## < Generating Alleles > ##########

    """
    Possible format.
    
    (1) Captioned or not (ex. "A*01:01:01:01" vs. "01:01:01:01")
    (2) with Double-Colon or not (ex. "A*01:01:01:01" vs. "A*01010101" / "106:01" vs. "10601")
    (3) Complete or not (ex. "A*01:01:01:01" vs. "A*01:01" / "01:01:01:01" vs. "01:01:01") 
    
    (4) as P-group or G-group ?
    
    """

    l_temp = []

    for i in range(0, len(HLA_names)):
    # for i in range(0, 1):

        curr_hla_name = HLA_names[i]
        curr_IAT_dict = IAT_dict[curr_hla_name]

        NumberOfAlleles = len(curr_IAT_dict)


        ### Picking alleles from *.iat table.

        RandomAlleles1 = [curr_IAT_dict.iat[random.randrange(0, NumberOfAlleles), random.choice(l_forChoice)] for j in range(0, _N)]
        RandomAlleles2 = [curr_IAT_dict.iat[random.randrange(0, NumberOfAlleles), random.choice(l_forChoice)] for j in range(0, _N)]

        print(std_MAIN_PROCESS_NAME + "Randomly generated alleles for {0}.\n".format(HLA_names[i]))
        print(RandomAlleles1)
        print(RandomAlleles2)


        ### Applying various options("Caption", "Incomplete", "NoDoubleColon")

        InCaseofGorPgroup = (_f_asPgroup or _f_asGgroup)

        RandomAlleles1 = pd.Series(RandomAlleles1).apply(lambda x : dressRandomFormat(curr_hla_name, x, InCaseofGorPgroup, _f_caption, _f_trimming, _f_DC))
        RandomAlleles2 = pd.Series(RandomAlleles2).apply(lambda x : dressRandomFormat(curr_hla_name, x, InCaseofGorPgroup, _f_caption, _f_trimming, _f_DC))


        print(std_MAIN_PROCESS_NAME + "Randomly generated alleles for {0}.\n".format(HLA_names[i]))
        print(RandomAlleles1)
        print(RandomAlleles2)


        ### Gathering each 2 colums to concatenate as DataFrame

        l_temp.append(RandomAlleles1)
        l_temp.append(RandomAlleles2)


    ### Generating ad DataFrame

    df_OUTPUT = pd.concat(l_temp, axis=1)
    print(df_OUTPUT.head())

    ### Making Dummy index

    l_temp = [['_'.join([header_ped[j], str(i)]) for i in range(0, _N)] for j in range(0, 4)]
    l_temp.append([random.randrange(0, 2) for i in range(0, _N)])
    l_temp.append([random.randrange(0, 2) for i in range(0, _N)])
    df_OUTPUT.index = pd.MultiIndex.from_arrays(l_temp, names=header_ped)


    ### Exporting OUTPUT

    print(std_MAIN_PROCESS_NAME + "Randomly generated *.ped file.\n")
    print(df_OUTPUT.head())

    # Export
    df_OUTPUT.to_csv(_out+".ped", sep='\t', header=False, index=True)


    return 0


def trimTheAllele(_allele):

    fields = _allele.split(':')

    if len(fields) < 3:
        return _allele
    elif len(fields) >= 3:

        [fields.pop() if random.choice([True, False]) else None for i in range(0, len(fields)-1)]
        return ':'.join(fields)



def removeDoubleColon(_allele):

    return re.sub(pattern=r'\:', repl="", string=_allele)


def dressRandomFormat(_hla_name, _given_allele,
                      _f_PorG, _f_caption, _f_trimming, _f_DC):

    ## Exception for No given allele
    if _given_allele == "0":
        # No need to dress given allele.
        return _given_allele

    #--------------#

    caption = ""
    Temp_Allele = _given_allele
    Trimmed_Allele = ""
    Trimmed_Allele_DC = ""

    ## Caption
    if _f_caption == 1:
        # Definitely add Gene caption.
        caption = _hla_name+"*"
    elif _f_caption == 2:
        caption = ""
    elif _f_caption == 3:
        # Randomly add Gene caption.
        caption = _hla_name+"*" if random.choice([True, False]) else ""

    if not _f_PorG:
        ## Trimming
        if _f_trimming == 1:
            # Definitely do trimming
            Trimmed_Allele = trimTheAllele(Temp_Allele)
        elif _f_trimming == 2:
            # No Trimming
            Trimmed_Allele = Temp_Allele
        elif _f_trimming == 3:
            # Randomly do trimming
            Trimmed_Allele = trimTheAllele(Temp_Allele) if random.choice([True, False]) else Temp_Allele
    else:
        # In case of P or G group, Skip trimming process.
        # (2018. 7. 1.) I think in case of "G-group" and "P-group", there would be no need to trimming.
        # (2018. 7. 2.) I missed some cases. For instance, "01:01:02" can appear even though we are dealing with P or G group.
        Trimmed_Allele = Temp_Allele

    ## No Double-Colon
    if _f_DC == 1:
        Trimmed_Allele_DC = Trimmed_Allele
    elif _f_DC == 2:
        # Definitely remove double-colon
        Trimmed_Allele_DC = removeDoubleColon(Trimmed_Allele)
    elif _f_DC == 3:
        # Randomly remove double-colon
        Trimmed_Allele_DC = removeDoubleColon(Trimmed_Allele) if random.choice([True, False]) else Trimmed_Allele



    return ''.join([caption, Trimmed_Allele_DC])





if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###############################################################################     
      
     < GeneratePED.py >
        
     Explanation.
        
        
     made by Wanson Choi.
     
    ###############################################################################
                                             '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-iat", help="\nInput \"*.iat\" file.\n\n", required=True)
    parser.add_argument("-o", help="\nOuput file prefix\n\n", required=True)
    parser.add_argument("-N", help="\nThe number of Samples(# of Rows)\n\n", required=True)

    # Flag for output format.
    output_format = parser.add_mutually_exclusive_group(required=True)
    output_format.add_argument("--as-Standard", help="\nOutput as only Standard 4-field alleles.\n\n", dest="onlyStandard", action='store_true')
    output_format.add_argument("--as-Ggroup", help="\nOutput as only G-group(No trimming off).\n\n", dest="onlyGgroup", action='store_true')
    output_format.add_argument("--as-Pgroup", help="\nOutput as only P-group(No trimming off).\n\n", dest="onlyPgroup", action='store_true')

    # Flag for additional features
    f_caption = parser.add_mutually_exclusive_group(required=False)
    f_caption.add_argument("--caption", help="\nCertainly add gene caption(ex. \"DRB1*\").\n\n", dest="Caption", action='store_true')
    f_caption.add_argument("--no-caption", help="\nNo gene caption(ex. \"DRB1*\").\n\n", dest="noCaption", action='store_true')
    # f_caption.add_argument("--random-caption", help="No gene caption(ex. \"DRB1*\").\n\n", dest="noCaption", action='store_true')

    f_trimming = parser.add_mutually_exclusive_group(required=False)
    f_trimming.add_argument("--trimming", help="\nCertainly do trimming off a few fields.\n\n", dest="Trim", action='store_true')
    f_trimming.add_argument("--no-trimming", help="\nNo trimming off a few fields.\n\n", dest="noTrim", action='store_true')
    # f_trimming.add_argument("--random-trimming", help="No trimming off a few fields.\n\n", dest="noTrim", action='store_true')

    f_DC = parser.add_mutually_exclusive_group(required=False)
    f_DC.add_argument("--doublecolon", help="\nCertainly add double-colons between fields.\n\n", dest="DC", action='store_true')
    f_DC.add_argument("--no-doublecolon", help="\nNo double-colons between fields.\n\n", dest="noDC", action='store_true')


    # <for Publish>
    args = parser.parse_args()

    # <for Testing>

    # # (Case 1) Classic combination of argument(All combination of "Standard", "G-group", "P-group", "Gene-Caption", "complete and incomplete", "Double-colon or not", etc.
    # args = parser.parse_args(["-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/DummyPED",
    #                           "-N", "30"])

    # # (Case 2) Only "Standard allele".
    # args = parser.parse_args(["-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/DummyPED",
    #                           "-N", "30",
    #                           "--as-Standard"])

    # # (Case 3) Only "Standard allele".
    # args = parser.parse_args(["-iat", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/DummyPED",
    #                           "-N", "200",
    #                           "--as-Standard",
    #                           "--no-caption",
    #                           "--doublecolon",
    #                           "--no-trimming"])

    print(args)


    ### Additional Flag processing

    """
    As three argument for "Caption", "Trimming" or "DoubleColon" are taken exclusively,
    I won't consider the case where both status are true.
    """

    ## Caption or Not.
    F_Caption = 0

    if args.Caption and not args.noCaption:
        F_Caption = 1 # add caption
    elif not args.Caption and args.noCaption:
        F_Caption = 2 # NO caption
    elif not args.Caption and not args.noCaption:
        F_Caption = 3 # Randomly

    ## Trimming or Not.
    F_Trmming = 0

    if args.Trim and not args.noTrim:
        F_Trmming = 1
    elif not args.Trim and args.noTrim:
        F_Trmming = 2
    elif not args.Trim and not args.noTrim:
        F_Trmming = 3

    ## Remove Double-colon or not.
    F_DC = 0

    if args.DC and not args.noDC:
        F_DC = 1
    elif not args.DC and args.noDC:
        F_DC = 2
    elif not args.DC and not args.noDC:
        F_DC = 3


    GeneratePED(args.iat, args.o, args.N,
                _f_asStandard=args.onlyStandard, _f_asGgroup=args.onlyGgroup, _f_asPgroup=args.onlyPgroup,
                _f_caption=F_Caption, _f_trimming=F_Trmming, _f_DC=F_DC)