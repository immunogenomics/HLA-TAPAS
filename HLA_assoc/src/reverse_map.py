# -*- coding: utf-8 -*-

import os, sys, re
from os.path import exists
import gzip

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
# p_HLA_marker = re.compile(r'^HLA_((\w+)\*[0-9A-Z:]+)$')
p_HLA_marker = re.compile(r'HLA_((\w+)\*\d{2,3}(:\d{2,3})+[A-Z]?)')


def reverse_map(_vcf, _out, _hped=None, _chped=None, _hat=None):

    """
    A module to perform reverse-mapping, from chped to hped, on a given bim file.
    (cf. NomenCleaner transforms hped to chped. e.g. G-group to 6-digit.)

    This is the function which Yang requested in the mail. (2020.04.26 00:22)
    """


    """
    Two major ways to get reverse-mapping information
    
    (1) hped and chped are given => directly make a mapping information from chped to hped(reverse way).
    (2) Just *.hat file is given => Basically same as performing NomenCleaner one more again.
    """

    ##### Preparing Reverse-mapping.

    __REVERSE_MAP__ = None

    if bool(_hped) and bool(_chped):

        ### (1) hped and chped are given.

        if not exists(_hped):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given HPED file('{}') can't be found.".format(_hped))
            sys.exit(-1)
        if not exists(_chped):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given CHPED file('{}') can't be found.".format(_chped))
            sys.exit(-1)


        __REVERSE_MAP__ = extract_reverse_map(_hped, _chped)


    elif bool(_hat):

        ### (2) Just *.hat file is given.

        # Not yet. (2020. 05. 01.)

        pass
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "None of HLA type mapping information given.")
        sys.exit(-1)



    ##### Main Reverse-mapping.

    flag_gzip = _vcf.endswith('.vcf.gz')

    with (gzip.open(_vcf, 'rt') if flag_gzip else open(_vcf, 'r')) as f_vcf, \
            gzip.open(_out, 'wb') as f_out, \
            open(_out + '.log', 'w') as f_log:

        count = 0

        for line in f_vcf:

            if line.startswith('##') or line.startswith('#'):
                # (1): Header lines
                # (2): Column line
                f_out.write(line.encode()) # Just write it out as it is.

            else:

                s = p_HLA_marker.search(line)

                if bool(s):
                    # (3-1): HLA markers (ex. HLA_A*01:01:01:01)

                    allele_label = s.group(1)   # A*01:01:01:01
                    hla = s.group(2)    # A

                    if allele_label in __REVERSE_MAP__[hla].keys():

                        allele_label_rev_mapped = __REVERSE_MAP__[hla][allele_label]
                        new_marker_label = 'HLA_{}'.format(allele_label_rev_mapped)
                        f_log.write("Before: {} => After: {}\n".format(s.group(0), new_marker_label))

                        new_line = p_HLA_marker.sub(new_marker_label, line, count=1)
                        f_out.write(new_line.encode())

                    else:
                        f_out.write(line.encode()) # Just write it out as it is.

                else:
                    # (3-2): the rest of (3-1). (ex. rs1234, HLA_A*01, etc.)
                    f_out.write(line.encode()) # Just write it out as it is.


            count += 1
            # if count > 2000 : break

    return _out



def extract_reverse_map(_hped, _chped):

    """
    Think hped as 'G-group' and chped as '6-digit'.

    """

    dict_HLA_reverse_map = {hla: {} for hla in HLA_names}

    with open(_hped, 'r') as f_hped, open(_chped, 'r') as f_chped:

        count = 0

        for line_hped, line_chped in zip(f_hped, f_chped):

            # print("\n{}".format(count))

            l_hped = re.split(r'\s+', line_hped.rstrip('\n'))
            # print("hped:\n{}".format(l_hped))

            l_chped = re.split(r'\s+', line_chped.rstrip('\n'))
            # print("chped:\n{}".format(l_chped))


            """
            0: FID
            1: IID
            2: PID
            3: MID
            4: Sex
            5: Phe
            
            6,7: HLA-A
            8,9: HLA-B
            ...
            20,21: HLA-DRB1
            
            """

            for i in range(len(HLA_names)):

                idx1 = 2*i + 6
                idx2 = idx1 + 1

                t_keys = dict_HLA_reverse_map[HLA_names[i]]

                if l_chped[idx1] not in t_keys:
                    dict_HLA_reverse_map[HLA_names[i]][l_chped[idx1]] = l_hped[idx1]

                if l_chped[idx2] not in t_keys:
                    dict_HLA_reverse_map[HLA_names[i]][l_chped[idx2]] = l_hped[idx2]



            count += 1
            # if count > 5 : break

        # Checking reverse-map dictionary
        # for k, v in dict_HLA_reverse_map.items():
        #     print("{}:\n{}".format(k, v))

    return dict_HLA_reverse_map



if __name__ == '__main__':

    [_bim, _out, _hped, _chped] = sys.argv[1:]
    reverse_map(_bim, _out, _hped, _chped)

