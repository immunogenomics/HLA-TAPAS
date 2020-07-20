# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
header_ped = ["FID", "IID", "PID", "MID", "Sex", "Phe"]

# Patterns
p_Suffix = re.compile(r'.+[A-Z]$')
p_Prefix = re.compile(r'^\w+\*')



class HATK_NomenCleaner(object):

    def __init__(self, _hped, _hat, _out, **kwargs):

        """

        """

        if not os.path.exists(_hped):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given HPED file('{}') doesn't exist.\n"
                                                "Please check the '--hped' argument again.\n"
                                                "".format(_hped))
            sys.exit()
        else:
            # Checking the number of columns
            f_hped = open(_hped, 'r')

            hped_1st_line = f_hped.readline()
            hped_1st_line = re.split(pattern='\s+', string=hped_1st_line.rstrip('\n'))

            if len(hped_1st_line) != 22:
                # # of columns of hped file must be 22(6 + 8*2).
                print(std_ERROR_MAIN_PROCESS_NAME + "The number of columns of given HPED file('{}') must be 22 but it is '{}'.\n"
                                                    "Please check the number of columns of that HPED file again."
                                                    "".format(_hped, len(hped_1st_line)))
                sys.exit()


        if not os.path.exists(_hat):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given HAT(HLA Allel Table) file('{}') doesn't exist.\n"
                                                "Please check the '--hat")



        self.chped = NomenCleaner(_hped, _hat, _out,
                                  __oneF=kwargs['__oneF'], __twoF=kwargs['__twoF'], __threeF=kwargs['__threeF'], __fourF=kwargs['__fourF'],
                                  __Ggroup=kwargs['__Ggroup'], __Pgroup=kwargs['__Pgroup'],
                                  __f_NoCaption=kwargs["__f_NoCaption"], __leave_NotFound=kwargs["__leave_NotFound"])




def NomenCleaner(_hped, _hat, _out, __f_NoCaption=False, __leave_NotFound=False, **kwargs):


    ##### < Output Format > #####

    FIELD = ''

    if kwargs['__oneF']:
        OUTPUT_FORMAT = 1
        FIELD = 'Maximum 1 field'
    elif kwargs['__twoF']:
        OUTPUT_FORMAT = 2
        FIELD = 'Maximum 2 fields'
    elif kwargs['__threeF']:
        OUTPUT_FORMAT = 3
        FIELD = 'Maximum 3 fields'
    elif kwargs['__fourF']:
        OUTPUT_FORMAT = 4
        FIELD = 'Maximum 4 fields'
    elif kwargs['__Ggroup']:
        OUTPUT_FORMAT = 5
        FIELD = 'G code'
    elif kwargs['__Pgroup']:
        OUTPUT_FORMAT = 6
        FIELD = 'P code'
    else:
        # Default mode
        OUTPUT_FORMAT = 0


    if OUTPUT_FORMAT != 0:
        print(std_MAIN_PROCESS_NAME + "Generating CHPED with {} HLA alleles.".format(FIELD))


    ##### < Loading Data > #####

    # *.hped
    __HPED__ = pd.read_csv(_hped, sep='\s+', header=None, dtype=str)
    # print("__HPED__ :\n{}\n".format(__HPED__.head()))

    # *.hat
    __HAT__ = pd.read_csv(_hat, sep='\s+', header=0, dtype=str, index_col=0)
    # print("__HAT__ :\n{}\n".format(__HAT__.head()))

    d__HAT__ = {HLA_names[i]: __HAT__.loc[HLA_names[i], :] for i in range(len(HLA_names))}

    # for k,v in d__HAT__.items():
    #     print("HLA : {}\n{}\n".format(k, v.head()))



    ##### < Main iteration > #####

    f_chped = open(_out+'.chped', 'w')
    f_chped_log = open(_out+'.chped.log', 'w')

    count = 0

    for eachRow in __HPED__.itertuples():

        # print(eachRow)

        [t_FID, t_IID] = eachRow[1:3]

        t_alleles = eachRow[7:]
        t_alleles = list(map(lambda x : p_Prefix.sub(string=x, repl=''), t_alleles)) # Ripping off gene prefix
        # print(t_alleles)

        ### [1] Allele conversion

        l_alleles_row = []

        for i in range(len(HLA_names)):

            idx1 = 2*i
            idx2 = 2*i + 1


            # chromosome 1
            if t_alleles[idx1] != '0':

                [t_converted_allele1, LOG_MESSAGE1] = getConvertedAllele2(HLA_names[i], t_alleles[idx1], d__HAT__[HLA_names[i]], OUTPUT_FORMAT, __leave_NotFound)

                # print("Converted Alleles : {}".format(t_converted_allele1))
                # print("LOG_MESSAGE :{}\n".format(LOG_MESSAGE1))

            else:

                t_converted_allele1 = t_alleles[idx1]
                LOG_MESSAGE1 = 'No conversion ({})'.format(t_alleles[idx1])




            # chromosome 2
            if t_alleles[idx2] != '0':

                [t_converted_allele2, LOG_MESSAGE2] = getConvertedAllele2(HLA_names[i], t_alleles[idx2], d__HAT__[HLA_names[i]], OUTPUT_FORMAT, __leave_NotFound)

                # print("Converted Alleles : {}".format(t_converted_allele2))
                # print("LOG_MESSAGE :{}\n".format(LOG_MESSAGE2))

            else:
                t_converted_allele2 = t_alleles[idx2]
                LOG_MESSAGE2 = 'No conversion ({})'.format(t_alleles[idx2])




            l_alleles_row.append(('%s*'%(HLA_names[i]) if (not __f_NoCaption and t_converted_allele1 != '0') else '') + t_converted_allele1)
            l_alleles_row.append(('%s*'%(HLA_names[i]) if (not __f_NoCaption and t_converted_allele2 != '0') else '') + t_converted_allele2)

            ### [2] logging the result
            f_chped_log.write('\t'.join([t_FID, t_IID, HLA_names[i]+'_1', LOG_MESSAGE1])+'\n')
            f_chped_log.write('\t'.join([t_FID, t_IID, HLA_names[i]+'_2', LOG_MESSAGE2])+'\n')



        ### Generating chped
        f_chped.write('\t'.join(list(eachRow[1:7]) + l_alleles_row) + '\n')


        count += 1
        # if count > 1 : break


    f_chped.close()
    f_chped_log.close()

    return _out+'.chped'





def getConvertedAllele2(_hla, _allele, _d__HAT__, _OUTPUT_FORMAT, __leave_NotFound):


    __RETURN__ = None # Converted allele
    LOG_MESSAGE = None


    ### < Setting `_from` > ###

    # `_from` := 'STANDARD', 'OLD', 'Ggroup' or 'Pgroup' column in HAT file.

    if _allele.endswith('G'):
        _from = 'Ggroup'
    elif _allele.endswith('P'):
        _from = 'Pgroup'
    else:
        if ':' in _allele:
            _from = 'STANDARD'
        else:
            _from = 'OLD'




    if ':' in _allele:

        # Colon (either STANDARD, G-group or P-group).

        [__RETURN__, LOG_MESSAGE] = get1stAllele2(_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)


    else:

        # No colon (either OLD or STANDARD).
        # Assumes it OLD first.

        if len(_allele) <= 3:

            # Exception handling for 1field allele (ex. DPB1*03 -> 0302(-> 102:01; OLD) / 03:01:01:01(STANDARD))
            # G-group and P-group alleles can't be here naturally. (ex. 0101P(01:01P) ... at least 5 characters.)

            flag_matched_STANDARD = _d__HAT__['STANDARD'].str.match(_allele)
            # Check both first.

            if flag_matched_STANDARD.any():

                # Give priority to 'STANDARD' over 'OLD'.

                [__RETURN__, LOG_MESSAGE] = get1stAllele2(_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)

            else:

                [__RETURN__, LOG_MESSAGE] = get1stAllele2(_allele, _d__HAT__, 'OLD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)




        else:

            ### Possibly OLD nomenclature

            # flag_matched_OLD = _d__HAT__[_from].str.match(_allele)

            flag_matched_OLD = _d__HAT__['OLD'].str.match(_allele).any()


            flag_matched_STANDARD = False


            ## Checking Suffix
            hasSuffix = p_Suffix.match(_allele)

            if hasSuffix:
                _allele2 = _allele[:-1]
                _allele_Suffix = _allele[-1]
            else:
                _allele2 = _allele
                _allele_Suffix = -1

            ### Main Digit-checking

            sr_TEMP = _d__HAT__['STANDARD'].str

            if len(_allele2) == 4:
                # 2 + 2 (+C) = 4

                trial1_allele = ':'.join([_allele2[:2], _allele2[2:4]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                Flag_trial1 = sr_TEMP.match(trial1_allele)

                if Flag_trial1.any():
                    [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    flag_matched_STANDARD = True
                else:
                    __RETURN__ = _allele if __leave_NotFound else "0"
                    if __leave_NotFound:
                        LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                    else:
                        LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

            elif len(_allele2) == 5:

                # (1) 2 + 3 (+C) = 5
                # (2) 3 + 2 (+C) = 5

                trial1_allele = ':'.join([_allele2[:2], _allele2[2:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                Flag_trial1 = sr_TEMP.match(trial1_allele)

                if Flag_trial1.any():
                    [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    flag_matched_STANDARD = True

                else:
                    trial2_allele = ':'.join([_allele2[:3], _allele2[3:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial2 = sr_TEMP.match(trial2_allele)

                    if Flag_trial2.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        flag_matched_STANDARD = True


                    else:
                        __RETURN__ = _allele if __leave_NotFound else "0"
                        if __leave_NotFound:
                            LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                        else:
                            LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

            elif len(_allele2) == 6:

                # (1) 2 + 2 + 2 (+C)
                # (2) 3 + 3 (+C)

                # (1) 2 + 2 + 2 (+C)
                trial1_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                Flag_trial1 = sr_TEMP.match(trial1_allele)

                if Flag_trial1.any():
                    [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    flag_matched_STANDARD = True
                else:
                    # (2) 3 + 3 (+C)
                    trial2_allele = ':'.join([_allele2[:3], _allele2[3:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial2 = sr_TEMP.match(trial2_allele)

                    if Flag_trial2.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        flag_matched_STANDARD = True
                    else:
                        __RETURN__ = _allele if __leave_NotFound else "0"
                        if __leave_NotFound:
                            LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                        else:
                            LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

            elif len(_allele2) == 7:

                # (1) 2 + 3 + 2 (+C) = 7
                # (2) 2 + 2 + 3 (+C) = 7
                # (3) 3 + 2 + 2 (+C) = 7

                # (1) 2 + 3 + 2 (+C) = 7
                trial1_allele = ':'.join([_allele2[:2], _allele2[2:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                Flag_trial1 = sr_TEMP.match(trial1_allele)

                if Flag_trial1.any():
                    [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    flag_matched_STANDARD = True

                else:
                    # (2) 2 + 2 + 3 (+C) = 7
                    trial2_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial2 = sr_TEMP.match(trial2_allele)

                    if Flag_trial2.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        flag_matched_STANDARD = True
                    else:
                        # (3) 3 + 2 + 2 (+C) = 7
                        trial3_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial3 = sr_TEMP.match(trial3_allele)

                        if Flag_trial3.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial3_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                            flag_matched_STANDARD = True
                        else:
                            __RETURN__ = _allele if __leave_NotFound else "0"
                            if __leave_NotFound:
                                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                            else:
                                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

            elif len(_allele2) == 8:

                # (1) 2 + 2 + 2 + 2 (+C)
                # (2) 3 + 2 + 3 (+C)
                # (3) 3 + 3 + 2 (+C)

                # (1) 2 + 2 + 2 + 2
                trial1_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:6], _allele2[6:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                Flag_trial1 = sr_TEMP.match(trial1_allele)

                if Flag_trial1.any():
                    [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    flag_matched_STANDARD = True
                else:
                    # (2) 3 + 2 + 3
                    trial2_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial2 = sr_TEMP.match(trial2_allele)

                    if Flag_trial2.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        flag_matched_STANDARD = True
                    else:
                        # (3) 3 + 3 + 2
                        trial3_allele = ':'.join([_allele2[:3], _allele2[3:6], _allele2[6:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial3 = sr_TEMP.match(trial3_allele)

                        if Flag_trial3.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial3_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                            flag_matched_STANDARD = True
                        else:
                            __RETURN__ = _allele if __leave_NotFound else "0"
                            if __leave_NotFound:
                                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                            else:
                                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

            elif len(_allele2) == 9:

                # (1) 2 + 3 + 2 + 2 (+C)
                # (2) 3 + 2 + 2 + 2 (+C)

                # (1) 2 + 3 + 2 + 2
                trial1_allele = ':'.join([_allele2[:2], _allele2[2:5], _allele2[5:7], _allele2[7:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                Flag_trial1 = sr_TEMP.match(trial1_allele)

                if Flag_trial1.any():
                    [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    flag_matched_STANDARD = True
                else:
                    # (2) 3 + 2 + 2 + 2
                    trial2_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:7], _allele2[7:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial2 = sr_TEMP.match(trial2_allele)

                    if Flag_trial2.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        flag_matched_STANDARD = True
                    else:
                        __RETURN__ = _allele if __leave_NotFound else "0"
                        if __leave_NotFound:
                            LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                        else:
                            LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

            elif len(_allele2) == 10:

                # (1) 3 + 3 + 2 + 2

                trial1_allele = ':'.join([_allele2[0:3], _allele2[3:6], _allele2[6:8], _allele2[8:10]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                Flag_trial1 = sr_TEMP.match(trial1_allele)

                if Flag_trial1.any():
                    [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, 'STANDARD', _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    flag_matched_STANDARD = True
                else:
                    __RETURN__ = _allele if __leave_NotFound else "0"
                    if __leave_NotFound:
                        LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                    else:
                        LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)





            if flag_matched_OLD and flag_matched_STANDARD:

                """
                (2020. 05. 15.)
                ex. HLA-B: 0702

                OLD -> STANDARD : 070202 -> 07:02:02    ... (1)
                STANDARD : 07:02:01:01                  ... (2) => This is proper.

                """

                # Need to check respective STANDARD names first.
                [return_old, log_message_old] = get1stAllele2(_allele, _d__HAT__, 'OLD', 4, __leave_NotFound, _allele0=_allele)
                [return_standard, log_message_standard] = [__RETURN__, LOG_MESSAGE]


                if return_old <= return_standard:

                    # (e.g.) 07:02:01:01 < 07:02:02

                    if _OUTPUT_FORMAT == 4:
                        return [return_old, log_message_old]
                    else:
                        _from = 'OLD'
                else:

                    if _OUTPUT_FORMAT == 4:
                        return [return_standard, log_message_standard]
                    else:
                        _from = 'STANDARD'

            elif flag_matched_OLD and not flag_matched_STANDARD:

                # Only found in OLD
                _from = 'OLD'

            elif not flag_matched_OLD and flag_matched_STANDARD:

                # Only found in STANDARD
                _from = 'STANDARD'

            else:

                # Wrong. Maybe not in the given IMGT version.
                if __leave_NotFound:
                    return [_allele, '[Leave Not Found] {} -> {}'.format(_allele, _allele)]
                else:
                    return ['0', "[Not Found] {} -> 0".format(_allele)]





            if _from == 'OLD':

                [__RETURN__, LOG_MESSAGE] = get1stAllele2(_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)




            else:

                ## Checking Suffix
                hasSuffix = p_Suffix.match(_allele)

                if hasSuffix:
                    _allele2 = _allele[:-1]
                    _allele_Suffix = _allele[-1]
                else:
                    _allele2 = _allele
                    _allele_Suffix = -1

                ### Main Digit-checking

                sr_TEMP = _d__HAT__[_from].str

                if len(_allele2) == 4:
                    # 2 + 2 (+C) = 4

                    trial1_allele = ':'.join([_allele2[:2], _allele2[2:4]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial1 = sr_TEMP.match(trial1_allele)

                    if Flag_trial1.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    else:
                        __RETURN__ = _allele if __leave_NotFound else "0"
                        if __leave_NotFound:
                            LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                        else:
                            LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                elif len(_allele2) == 5:

                    # (1) 2 + 3 (+C) = 5
                    # (2) 3 + 2 (+C) = 5

                    trial1_allele = ':'.join([_allele2[:2], _allele2[2:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial1 = sr_TEMP.match(trial1_allele)

                    if Flag_trial1.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    else:
                        trial2_allele = ':'.join([_allele2[:3], _allele2[3:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial2 = sr_TEMP.match(trial2_allele)

                        if Flag_trial2.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        else:
                            __RETURN__ = _allele if __leave_NotFound else "0"
                            if __leave_NotFound:
                                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                            else:
                                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                elif len(_allele2) == 6:

                    # (1) 2 + 2 + 2 (+C)
                    # (2) 3 + 3 (+C)

                    # (1) 2 + 2 + 2 (+C)
                    trial1_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial1 = sr_TEMP.match(trial1_allele)

                    if Flag_trial1.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    else:
                        # (2) 3 + 3 (+C)
                        trial2_allele = ':'.join([_allele2[:3], _allele2[3:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial2 = sr_TEMP.match(trial2_allele)

                        if Flag_trial2.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        else:
                            __RETURN__ = _allele if __leave_NotFound else "0"
                            if __leave_NotFound:
                                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                            else:
                                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                elif len(_allele2) == 7:

                    # (1) 2 + 3 + 2 (+C) = 7
                    # (2) 2 + 2 + 3 (+C) = 7
                    # (3) 3 + 2 + 2 (+C) = 7

                    # (1) 2 + 3 + 2 (+C) = 7
                    trial1_allele = ':'.join([_allele2[:2], _allele2[2:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial1 = sr_TEMP.match(trial1_allele)

                    if Flag_trial1.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    else:
                        # (2) 2 + 2 + 3 (+C) = 7
                        trial2_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial2 = sr_TEMP.match(trial2_allele)

                        if Flag_trial2.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        else:
                            # (3) 3 + 2 + 2 (+C) = 7
                            trial3_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                            Flag_trial3 = sr_TEMP.match(trial3_allele)

                            if Flag_trial3.any():
                                [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial3_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                            else:
                                __RETURN__ = _allele if __leave_NotFound else "0"
                                if __leave_NotFound:
                                    LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                                else:
                                    LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                elif len(_allele2) == 8:

                    # (1) 2 + 2 + 2 + 2 (+C)
                    # (2) 3 + 2 + 3 (+C)
                    # (3) 3 + 3 + 2 (+C)

                    # (1) 2 + 2 + 2 + 2
                    trial1_allele = ':'.join([_allele2[:2], _allele2[2:4], _allele2[4:6], _allele2[6:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial1 = sr_TEMP.match(trial1_allele)

                    if Flag_trial1.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    else:
                        # (2) 3 + 2 + 3
                        trial2_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial2 = sr_TEMP.match(trial2_allele)

                        if Flag_trial2.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        else:
                            # (3) 3 + 3 + 2
                            trial3_allele = ':'.join([_allele2[:3], _allele2[3:6], _allele2[6:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                            Flag_trial3 = sr_TEMP.match(trial3_allele)

                            if Flag_trial3.any():
                                [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial3_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                            else:
                                __RETURN__ = _allele if __leave_NotFound else "0"
                                if __leave_NotFound:
                                    LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                                else:
                                    LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                elif len(_allele2) == 9:

                    # (1) 2 + 3 + 2 + 2 (+C)
                    # (2) 3 + 2 + 2 + 2 (+C)

                    # (1) 2 + 3 + 2 + 2
                    trial1_allele = ':'.join([_allele2[:2], _allele2[2:5], _allele2[5:7], _allele2[7:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial1 = sr_TEMP.match(trial1_allele)

                    if Flag_trial1.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    else:
                        # (2) 3 + 2 + 2 + 2
                        trial2_allele = ':'.join([_allele2[:3], _allele2[3:5], _allele2[5:7], _allele2[7:]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                        Flag_trial2 = sr_TEMP.match(trial2_allele)

                        if Flag_trial2.any():
                            [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial2_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                        else:
                            __RETURN__ = _allele if __leave_NotFound else "0"
                            if __leave_NotFound:
                                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                            else:
                                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)


                elif len(_allele2) == 10:

                    # (1) 3 + 3 + 2 + 2

                    trial1_allele = ':'.join([_allele2[0:3], _allele2[3:6], _allele2[6:8], _allele2[8:10]]) + (_allele_Suffix if _allele_Suffix != -1 else "")
                    Flag_trial1 = sr_TEMP.match(trial1_allele)

                    if Flag_trial1.any():
                        [__RETURN__, LOG_MESSAGE] = get1stAllele2(trial1_allele, _d__HAT__, _from, _OUTPUT_FORMAT, __leave_NotFound, _allele0=_allele)
                    else:
                        __RETURN__ = _allele if __leave_NotFound else "0"
                        if __leave_NotFound:
                            LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
                        else:
                            LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)




    return [__RETURN__, LOG_MESSAGE]




def get1stAllele2(_allele, _df, _from, _OUTPUT_FORMAT, __leave_NotFound=False, _allele0=None):


    __1st_Allele__ = None
    LOG_MESSAGE = None

    sr_from = _df[_from]


    ##### < HAT column to use (`_to`) > #####

    if _OUTPUT_FORMAT == 0:

        # Default mode. (`_from` == `_to`)

        if _from == 'STANDARD':
            _to = 1
        elif _from == 'OLD':
            _to = 1
        elif _from == 'Ggroup':
            _to = 3
        elif _from == 'Pgroup':
            _to = 4
        else:
            _to = -1

    elif 0 < _OUTPUT_FORMAT <= 4:
        _to = 1
    elif _OUTPUT_FORMAT == 5:
        _to = 3
    elif _OUTPUT_FORMAT == 6:
        _to = 4
    else:
        _to = -1



    ##### < Main Matching job > #####

    ## (1) Exact match
    Flag_exact_match = sr_from.str.match(_allele + '$')
    f_exact = False

    if Flag_exact_match.any():

        # (ex. '03:28' -> '03:28')

        df_candidate = _df.loc[Flag_exact_match, :]
        __1st_Allele__ = df_candidate.iat[0, _to]

        f_exact = True

    else:

        ## (2) Just match
        Flag_match = sr_from.str.match(_allele)

        if Flag_match.any():

            # (ex. '03:28' -> '03:280')

            df_candidate = _df.loc[Flag_match, :]
            __1st_Allele__ = df_candidate.iat[0, _to]

        else:

            ## (3) Nothing matched.
            __1st_Allele__ = _allele if __leave_NotFound else '0'

            if __leave_NotFound:
                LOG_MESSAGE = '[Leave Not Found] {} -> {}'.format(_allele, _allele)
            else:
                LOG_MESSAGE = "[Not Found] {} -> 0".format(_allele)

            return [__1st_Allele__, LOG_MESSAGE]




    ##### < Log message & Field trimming > #####


    if _from == 'OLD':

        if _to == 1:

            if _OUTPUT_FORMAT == 4:

                # to Standard
                LOG_MESSAGE = "[{}] {} -> {}(OLD) => {}(STANDARD)".format(
                    'Exact match' if f_exact else 'Approximated', _allele, df_candidate.iat[0, 2], __1st_Allele__)

                if not f_exact and df_candidate.shape[0] > 0:
                    LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                        df_candidate.iloc[:, [2, 1]].apply(' -> '.join, axis=1).tolist())])

            else:
                # 0 <= _OUTPUT_FORMAT < 4

                __1st_Allele_cutted__ = FieldCutter(':'.join(a+b for a,b in zip(_allele0[::2], _allele0[1::2])), __1st_Allele__, _OUTPUT_FORMAT)

                # to Standard
                LOG_MESSAGE = "[{}] {} -> ({}(OLD) -> {}(STANDARD)) => {}".format(
                    'Exact match' if f_exact else 'Approximated', _allele, df_candidate.iat[0, 2], __1st_Allele__, __1st_Allele_cutted__)

                if not f_exact and df_candidate.shape[0] > 0:
                    LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                        df_candidate.iloc[:, [2, 1]].apply(' -> '.join, axis=1).tolist())])


                __1st_Allele__ = __1st_Allele_cutted__



        elif _to == 3 or _to == 4:
            # to Ggroup or Pgroup
            LOG_MESSAGE = "[{}] {} -> ({}(OLD) -> {}(STANDARD)) => {}".format(
                'Exact match' if f_exact else 'Approximated', _allele, df_candidate.iat[0, 2], df_candidate.iat[0, 1], __1st_Allele__)

            if not f_exact and df_candidate.shape[0] > 0:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                    df_candidate.iloc[:, [2, 1, _to]].apply(' -> '.join, axis=1).tolist())])



    elif _from == 'Ggroup':

        s_temp = '({}'.format(_allele) if _allele0 == _allele else '{} -> ({}'.format(_allele0, _allele)

        if _to == 1:

            # There would be no duplicates in this case block.
            f_exact_Ggroup = f_exact and (df_candidate.shape[0] == 1)

            if _OUTPUT_FORMAT == 4:

                # to Standard
                LOG_MESSAGE = "[{}] {} => {}(STANDARD))".format(
                    'Exact match' if f_exact_Ggroup else 'Approximated', s_temp, __1st_Allele__)

                if not f_exact_Ggroup:
                    LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, _to].tolist())])

            else:
                # 0 < _OUTPUT_FORMAT < 4
                # the case "_OUTPUT_FORMAT == 0" can't be here. It will go to next "_to == 3" case block.

                __1st_Allele_cutted__ = FieldCutter(_allele, __1st_Allele__, _OUTPUT_FORMAT)

                LOG_MESSAGE = "[{}] {} -> {}(STANDARD)) => {}".format(
                    'Exact match' if f_exact_Ggroup else 'Approximated', s_temp, __1st_Allele__, __1st_Allele_cutted__)

                if not f_exact_Ggroup:
                    LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, _to].tolist())])


                __1st_Allele__ = __1st_Allele_cutted__



        elif _to == 3:

            # to G-group itself.

            df_candidate = df_candidate.drop_duplicates(['Ggroup'])
            # print(df_candidate)
            f_exact_Ggroup = f_exact and (df_candidate.shape[0] == 1)

            LOG_MESSAGE = "[{}] ({} => {})".format(
                'Exact match' if f_exact_Ggroup else 'Approximated', _allele0, __1st_Allele__)

            if not f_exact_Ggroup:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, _to].tolist())])



        elif _to == 4:
            # to Pgroup

            df_candidate = df_candidate.drop_duplicates(['Ggroup', 'Pgroup'])
            # print(df_candidate)
            f_exact_Ggroup = f_exact and (df_candidate.shape[0] == 1)

            LOG_MESSAGE = "[{}] {} => {})".format(
                'Exact match' if f_exact else 'Approximated', s_temp, __1st_Allele__)

            if not f_exact_Ggroup:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                    df_candidate.iloc[:, [3, _to]].apply(' -> '.join, axis=1).tolist())])



    elif _from == 'Pgroup':

        s_temp = '({}'.format(_allele) if _allele0 == _allele else '{} -> ({}'.format(_allele0, _allele)


        if _to == 1:

            f_exact_Pgroup = f_exact and (df_candidate.shape[0] == 1)

            if _OUTPUT_FORMAT == 4:

                # to Standard.
                LOG_MESSAGE = "[{}] {} => {}(STANDARD))".format(
                    'Exact match' if f_exact_Pgroup else 'Approximated', s_temp, __1st_Allele__)

                if not f_exact_Pgroup:
                    LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, _to].tolist())])

            else:
                # 0 < _OUTPUT_FORMAT < 4
                # the case "_OUTPUT_FORMAT == 0" can't be here. It will go to next "_to == 4" case block.

                __1st_Allele_cutted__ = FieldCutter(_allele, __1st_Allele__, _OUTPUT_FORMAT)

                LOG_MESSAGE = "[{}] {} -> {}(STANDARD)) => {}".format(
                    'Exact match' if f_exact_Pgroup else 'Approximated', s_temp, __1st_Allele__, __1st_Allele_cutted__)

                if not f_exact_Pgroup:
                    LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, _to].tolist())])


                __1st_Allele__ = __1st_Allele_cutted__



        elif _to == 3:
            # to G-group

            df_candidate = df_candidate.drop_duplicates(['Ggroup', 'Pgroup'])
            f_exact_Pgroup = f_exact and (df_candidate.shape[0] == 1)

            LOG_MESSAGE = "[{}] {} => {})".format(
                'Exact match' if f_exact_Pgroup else 'Approximated', s_temp, __1st_Allele__)

            if not f_exact_Pgroup:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                    df_candidate.iloc[:, [4, _to]].apply(' -> '.join, axis=1).tolist())])



        elif _to == 4:
            # to P-group itself.

            df_candidate = df_candidate.drop_duplicates(['Pgroup'])
            f_exact_Pgroup = f_exact and (df_candidate.shape[0] == 1)

            LOG_MESSAGE = "[{}] ({} => {})".format(
                'Exact match' if f_exact_Pgroup else 'Approximated', _allele0, __1st_Allele__)

            if not f_exact_Pgroup:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, _to].tolist())])



    else:

        if _to == 1:

            if _OUTPUT_FORMAT == 4:

                # to Standard
                LOG_MESSAGE = "[{}] {} => {}(STANDARD)".format(
                    'Exact match' if f_exact else 'Approximated', _allele0 if bool(_allele0) else _allele, __1st_Allele__)

                if not f_exact and df_candidate.shape[0] > 0:
                    LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, 1].tolist())])

            else:
                # 0 <= _OUTPUT_FORMAT < 4

                __1st_Allele_cutted__ = FieldCutter(_allele, __1st_Allele__, _OUTPUT_FORMAT)

                # to Standard
                LOG_MESSAGE = "[{}] ({} -> {}(STANDARD)) => {}".format(
                    'Exact match' if f_exact else 'Approximated', _allele0 if bool(_allele0) else _allele, __1st_Allele__, __1st_Allele_cutted__)

                if not f_exact and df_candidate.shape[0] > 0:
                    LOG_MESSAGE = ' '.join(
                        [LOG_MESSAGE, "(Possible candidates: {})".format(df_candidate.iloc[:, 1].tolist())])


                __1st_Allele__ = __1st_Allele_cutted__



        elif _to == 3 or _to == 4:
            # to Ggroup
            LOG_MESSAGE = "[{}] {} -> ({} => {})".format(
                'Exact match' if f_exact else 'Approximated', _allele0 if bool(_allele0) else _allele, df_candidate.iat[0, 1], __1st_Allele__)

            if not f_exact and df_candidate.shape[0] > 0:
                LOG_MESSAGE = ' '.join([LOG_MESSAGE, "(Possible candidates: {})".format(
                    df_candidate.iloc[:, [1, _to]].apply(' -> '.join, axis=1).tolist())])




    return [__1st_Allele__, LOG_MESSAGE]





def FieldCutter(_allele0, _1st_allele, _field_format):

    if _field_format > 4:
        return _allele0


    if 0 < _field_format <= 4:
        N_field = _field_format
    else:
        N_field = len(_allele0.split(':'))


    l_allele_found = _1st_allele.split(':')


    while N_field < len(l_allele_found) and len(l_allele_found) > 0:
        l_allele_found.pop()


    __RETURN__ = ':'.join(l_allele_found)
    # print("OUTPUT : {}".format(__RETURN__))


    return __RETURN__





# if __name__ == '__main__':
#
#
#     parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
#                                      description=textwrap.dedent('''\
#     #########################################################################################
#
#         < NomenCleaner.py >
#
#         - Transforms *.hped file to *.chped file.
#         - *.hat file must be given as input file.
#
#
#
#
#
#
#     #########################################################################################
#                                      '''),
#                                      add_help=False)
#
#     parser._optionals.title = "OPTIONS"
#
#     parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
#
#     # Input (1) : *.ped file
#     PED_TYPE = parser.add_mutually_exclusive_group(required=True)
#     PED_TYPE.add_argument("--hped", help="\nHLA Type Data with raw HLA allele(ex. 0101).\n\n", dest="hped")
#
#     # Input (2) : *.iat file
#     parser.add_argument("-hat", help="\nHLA Allele Table file(*.hat).\n\n", required=True)
#     parser.add_argument("-imgt", help="\nSpecifying the IMGT-HLA version.\n\n", required=True)
#
#     # Ouptut Prefix
#     parser.add_argument("--out", "-o", help="\nOutput file prefix.\n\n", required=True)
#     parser.add_argument("--leave-NotFound",
#                         help="\nLeaving HLA alleles which can't be found in given *.hat file(Novel or Erroneous allele) intact.\n\n",
#                         action='store_true')
#
#     # Output format
#     output_digit_selection = parser.add_mutually_exclusive_group()
#     output_digit_selection.add_argument("--1field", help="\nMake converted HLA alleles have maximum 1 field.\n\n",
#                                         action="store_true", dest="oneF")
#     output_digit_selection.add_argument("--2field", help="\nMake converted HLA alleles have maximum 2 fields.\n\n",
#                                         action="store_true", dest="twoF")
#     output_digit_selection.add_argument("--3field", help="\nMake converted HLA alleles have maximum 3 fields.\n\n",
#                                         action="store_true", dest="threeF")
#     output_digit_selection.add_argument("--4field", help="\nMake converted HLA alleles have maximum 4 fields\n\n",
#                                         action="store_true", dest="fourF")
#     output_digit_selection.add_argument("--Ggroup", help="\nMake converted HLA alleles have G code names.\n\n",
#                                         action="store_true")
#     output_digit_selection.add_argument("--Pgroup", help="\nMake converted HLA alleles have P code names.\n\n",
#                                         action="store_true")
#
#     # Flag to remove HLA gene caption.
#     parser.add_argument("--NoCaption", help="\nMake converted HLA alleles NOT have HLA gene prefix(ex. \"A*\").\n\n", action='store_true')
#
#
#
#     ##### <for Test> #####
#
#     # in OS X
#     # _hped = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.ruined.nocolon.hped'
#     # _hat = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/HLA_ALLELE_TABLE.imgt3320.hat'
#     # _out = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.ruined.nocolon.20190927'
#
#     # in Ubuntu
#     # _hped = '/home/wanson/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.ruined.hped'
#     # _hat = '/home/wanson/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/HLA_ALLELE_TABLE.imgt3320.hat'
#     # _out = '/home/wanson/Dropbox/_Sync_MyLaptop/Projects/HATK/NomenCleanerv3/dummy_freeze/DummyCHPED.10.ruined.20190922'
#
#     # for _output_format in ['DEFAULT', '1field', '2field', '3field', '4field', 'Ggroup', 'Pgroup']:
#     #
#     #     OUT = _out + '.{}'.format(_output_format.upper())
#     #
#     #     if _output_format == 'DEFAULT':
#     #         args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', OUT, '-imgt', '3320', "--NoCaption"]) # No field format given.
#     #     else:
#     #         args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', OUT, '-imgt', '3320', "--NoCaption", '--{}'.format(_output_format)]) # No field format given.
#     #
#     #     # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out, '-imgt', '3320', "--NoCaption", "--4field"]) # No field format given.
#     #     # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out, '-imgt', '3320', "--NoCaption", "--3field"]) # No field format given.
#     #     # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out, '-imgt', '3320', "--NoCaption", "--2field"]) # No field format given.
#     #     # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out, '-imgt', '3320', "--NoCaption", "--1field"]) # No field format given.
#     #     # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out, '-imgt', '3320', "--NoCaption", "--Ggroup"]) # No field format given.
#     #     # args = parser.parse_args(['--hped', _hped, '-hat', _hat, '-o', _out, '-imgt', '3320', "--NoCaption", "--Pgroup"]) # No field format given.
#     #
#     #     ##### <for Publication> #####
#     #
#     #     # args = parser.parse_args()
#     #     print(args)
#     #
#     #
#     #     NomenCleaner(args.hped, args.hat, args.imgt, args.out,
#     #                  __f_NoCaption=args.NoCaption, __leave_NotFound=args.leave_NotFound,
#     #                  __oneF=args.oneF, __twoF=args.twoF, __threeF=args.threeF, __fourF=args.fourF, __Ggroup=args.Ggroup, __Pgroup=args.Pgroup)
#
#
#
#     ##### <for Publication> #####
#
#     args = parser.parse_args()
#     print(args)
#
#
#     NomenCleaner(args.hped, args.hat, args.imgt, args.out,
#                  __f_NoCaption=args.NoCaption, __leave_NotFound=args.leave_NotFound,
#                  __oneF=args.oneF, __twoF=args.twoF, __threeF=args.threeF, __fourF=args.fourF, __Ggroup=args.Ggroup, __Pgroup=args.Pgroup)
