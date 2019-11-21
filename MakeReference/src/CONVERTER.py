#!/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

- Created by Wansun Choi, 2017/12/06
- CONVERTER.py

SPECIFICATION : This program converts results files from HLA Imputation software 'HIBAG' or 'Axiom HLA Typing Analysis Services' to .ped file that can be used in SNP2HLA

INPUT file : Result Files from 'HIBAG(https://bioconductor.org/packages/release/bioc/html/HIBAG.html)' or 'Axiom HLA Typing Analysis Services(https://www.thermofisher.com/order/catalog/product/000.911)'
OUTPUT file : .ped file

USAGE : ./COOKFORMATS.py PLATFORM(HIBAG/AXIOM) OUTPUT_PREFIX RESULTFILES(Multiple Files or one Folder)


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


import os, sys, re
import pandas as pd
import argparse, textwrap
from pathlib import Path

# main() function
def CONVERTER(_PLATFORM, _INPUT, _OUTPUT = "CONVERTED"):

    print("{0} {1} {2}".format(_PLATFORM, _INPUT, _OUTPUT))

    ########## <Core Variables> ##########

    Broad_Order = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
    # HLA_names = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
    # Broad_Order = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

    filelist = []

    HLA_patterns = {'A': re.compile('\.A|_A_|HLA-A|HLA_A'),
                    'C': re.compile('\.C|_C_|HLA-C|HLA_C'),
                    'B': re.compile('\.B|_B_|HLA-B|HLA_B'),
                    'DRB1': re.compile('\.DRB1|_DRB1_|HLA-DRB1|HLA_DRB1'),
                    'DQA1': re.compile('\.DQA1|_DQA1_|HLA-DQA1|HLA_DQA1'),
                    'DQB1': re.compile('\.DQB1|_DQB1_|HLA-DQB1|HLA_DQB1'),
                    'DPA1': re.compile('\.DPA1|_DPA1_|HLA-DPA1|HLA_DPA1'),
                    'DPB1': re.compile('\.DPB1|_DPB1_|HLA-DPB1|HLA_DPB1')
                    }

    HLA_FILENAMES = {'A': None,
                     'C': None,
                     'B': None,
                     'DRB1': None,
                     'DQA1': None,
                     'DQB1': None,
                     'DPA1': None,
                     'DPB1': None
                     }

    HLA_DataFrame = {'A': pd.DataFrame(),
                     'C': pd.DataFrame(),
                     'B': pd.DataFrame(),
                     'DRB1': pd.DataFrame(),
                     'DQA1': pd.DataFrame(),
                     'DQB1': pd.DataFrame(),
                     'DPA1': pd.DataFrame(),
                     'DPB1': pd.DataFrame()
                     }

    isDIR = 0
    filelist = pd.Series()
    intermediate_dir = pd.Series()

    FIRST_APPEAR = -1



    ########## <Control Flags> ##########

    CHEKCING_INPUT_TYPE = 1
    CHEKCING_FILENAMES = 1
    CONVERTING = 1


    ########## <INPUT type case classification> ##########

    if CHEKCING_INPUT_TYPE:

        """
        (1) 폴더로 주거나,
        (2) multiple files로 주거나.
        
        이걸 잘 나눠야 할거같은데, 우선 폴더로 주는 경우 _INPUT이 단일 string이거나 length가 1인 리스트라고 받고,
        multiple file의 경우 반드시 리스트라고 가정할 것.
        
        귀찮다. element의 개수가 1개이건 그 이상이건 일단 무조건 list형태로 받는가고 가정하자.
        
        
        (2018.4.3) 예를 들어서 그냥 multiple files로 받는다고 했을때, 이게 꼭 HLA_names순으로 받을거라는 보장이 없잖슴. 그렇기 때문에
        intermediate_dir과 filelist는 HLA_name순서하고 상관 없이 얘네 둘 사이의 순서만 매칭되는 거임. 
        
        """

        if (type(_INPUT) == list) or (type(_INPUT) == tuple):

            if os.path.isdir(_INPUT[-1]):
                print("\n\nGiven _INPUT is a dir.\n")

                isDIR = 1
                intermediate_dir = pd.Series(_INPUT[-1])
                filelist = pd.Series(os.listdir(_INPUT[-1]))

            else:
                # multiple file들인 경우
                print("\n\nGiven _INPUT are files\n")

                filelist = pd.Series(_INPUT)



        elif type(_INPUT) == str:
            print("\n\nGiven _INPUT is a single string\n")
            # single file or directory


        print("\n\nreal file name(s) :\n{0}".format(filelist))
        print("\n\nintermediate_path :\n{0}".format(intermediate_dir))



    if CHEKCING_FILENAMES:

        ########## <Checking Filenames> ##########

        for i in range(0, len(Broad_Order)):

            # 이름이 나타나는 boolean series
            if isDIR:
                sr_bool = filelist.apply(lambda x : True if re.search(pattern=HLA_patterns[Broad_Order[i]], string=x) else False)
            else:
                sr_bool = filelist.apply(lambda x: Path(x).name).apply(lambda x : True if re.search(pattern=HLA_patterns[Broad_Order[i]], string=x) else False)

            sr_found = filelist[sr_bool]

            if len(sr_found) > 0:
                # 찾은게 있는 경우
                # print(sr_found)

                HLA_FILENAMES[Broad_Order[i]] = os.path.join(intermediate_dir.iat[0], sr_found.iat[0]) if isDIR else sr_found.iat[0]


        print("\n\n<Found each of HLA files>")
        for k, v in HLA_FILENAMES.items():
            print("{0} : {1}".format(k, v))


    if CONVERTING:

        ########## <Main Job> ##########

        if _PLATFORM == "AXIOM":

            print("\n\nPerforming Converting job for AXIOM.\n")

            ### Loading files as DataFrames

            for i in range(0, len(Broad_Order)):

                if HLA_FILENAMES[Broad_Order[i]]: # 파일 이름을 잘 찾은 경우
                    HLA_DataFrame[Broad_Order[i]] = pd.read_table(filepath_or_buffer=HLA_FILENAMES[Broad_Order[i]], sep=' |\t',
                                                                engine='python', header=None, index_col=0, dtype=str, names=['Sample_ID', 'pair', 'Labels', 'p1', 'p2'])

                    if HLA_DataFrame[Broad_Order[i]].shape[0] % 2:
                        # Axiom의 경우 각 HLA 파일별로 row의 개수가 반드시 짝수개여야함. 혹시 몰라서 예외처리 해놓음
                        print("File for HLA_{0} has odd number of row. Please check again!".format(Broad_Order[i]))
                        sys.exit()

                else:
                    print("file for {0} doesn't exist!".format(Broad_Order[i]))


            # # 로드된 데이터프레임 결과물 확인
            # for i in range(0, len(HLA_names)):
            #     print("\nHLA_{0}".format(HLA_names[i]))
            #     print(HLA_DataFrame[HLA_names[i]].head())



            # 첫 번째로 나타나는 HLA index 구하기
            for i in range(0, len(Broad_Order)):
                if not HLA_DataFrame[Broad_Order[i]].empty:
                    FIRST_APPEAR = i
                    break # 이렇게 break을 했으니 가장 첫 번째 나타나는 데이터프레임

                # 위에서 row수 짝수인건 확인했으니 여기선 스킵하자


            """
            사실 label가져오는 작업은 굉장히 중요하고, 사실 Series로 만들어서 unique함수 써버릴까 생각했는데 이 상황에선 오히려 이렇게 하는게 낫겠음.
            똑바로 주어진다는 가정하에 저렇게 unique를 시키든 짝수,홀수번만 캐오든 어쨌든 아무 문제 없어야함. 데이터가 딱 저 상태인게 맞다면
            """


            if not (FIRST_APPEAR < 0):

                """
                가끔씩 HLA 한 두개 삔꾸난다는 가정하에, 이런 애들은 HLA_DataFrame에서 empty DataFrame으로 존재할텐데 여기다가 바로 pivot함수쓰면 에러나서 안될듯.
                그냥 for문을 한번 돌면서 처리하는게 나을거같음.
                """

                # pivot함수 만들기 앞서 Sample_id 집합 구하기

                ped_Sample_id = HLA_DataFrame[Broad_Order[FIRST_APPEAR]].pivot(columns='pair', values='Labels').index.tolist()


                l_AXIOM = [] # HLA별로 처리하고 나온 Series의 리스트

                ### 각 HLA별 allele name따오기

                for i in range(0, len(Broad_Order)):
                    if not HLA_DataFrame[Broad_Order[i]].empty:
                        # 데이터프레임이 정상적으로 데이터가 있어서 읽어들여놓은 경우 => pivot함수 맘놓고 사용
                        l_AXIOM.append(HLA_DataFrame[Broad_Order[i]].pivot(columns='pair', values='Labels'))

                    else:
                        # 0으로 2컬럼 만들어서 채워줘야함.
                        l_NULL = ['0' for i in range(0, len(ped_Sample_id))]
                        sr_NULL = pd.concat([pd.Series(l_NULL), pd.Series(l_NULL)], axis=1)
                        sr_NULL.index = pd.Index(ped_Sample_id)

                        l_AXIOM.append(sr_NULL)


                ped_Family_ID = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Paternal_ID = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Maternal_ID = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Sex = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Affection = ['0' for i in range(0, len(ped_Sample_id))]

                #### *** 짱중요, 여기서 Sex, Affection어떻게 반영할지 물어볼것(이거 파싱하는 역할만 따로 추가해야할지 모름)


                df_AXIOM = pd.concat(l_AXIOM, axis=1)
                df_AXIOM.index = pd.MultiIndex.from_arrays([ped_Family_ID, ped_Sample_id, ped_Paternal_ID, ped_Maternal_ID, ped_Sex, ped_Affection])

                print(df_AXIOM.head())

                df_AXIOM.to_csv(_OUTPUT+'.axiom.ped', sep='\t', header=False, index=True)


            else:
                # 처음 초기화 해놓은 -1 그대로 가지고 있는 상황 <=> HLA어느 파일도 못찾은 상황
                print("No Axiom file for any HLA! Please check input file directory and files again.")
                sys.exit()




        elif _PLATFORM == "HIBAG":

            print("\n\nPerforming Converting job for HIBAG.\n")

            ### Loading files as DataFrames

            for i in range(0, len(Broad_Order)):

                if HLA_FILENAMES[Broad_Order[i]]: # 파일 이름을 잘 찾은 경우
                    HLA_DataFrame[Broad_Order[i]] = pd.read_table(filepath_or_buffer=HLA_FILENAMES[Broad_Order[i]], sep=' |\t',
                                                                engine='python', header=None, index_col=0, dtype=str,
                                                                names=['Sample_ID', 'al1', 'al2', 'p'],
                                                                usecols=[0,1,2])

                else:
                    print("file for {0} doesn't exist!".format(Broad_Order[i]))



            # 첫 번째로 나타나는 HLA index 구하기
            for i in range(0, len(Broad_Order)):
                if not HLA_DataFrame[Broad_Order[i]].empty:
                    FIRST_APPEAR = i
                    break # 이렇게 break을 했으니 가장 첫 번째 나타나는 데이터프레임



            if not (FIRST_APPEAR < 0):
                # 정상적으로 적어도 하나는 찾은 경우

                ped_Sample_id = HLA_DataFrame[Broad_Order[FIRST_APPEAR]].index.tolist()

                l_HIBAG = []

                for i in range(0, len(Broad_Order)):

                    if not HLA_DataFrame[Broad_Order[i]].empty:

                        print("DataFrame for HLA_{0}:\n{1}".format(Broad_Order[i], HLA_DataFrame[Broad_Order[i]].head()))

                        l_HIBAG.append(HLA_DataFrame[Broad_Order[i]].applymap(lambda x : re.sub(pattern='\:', repl='', string=x)))

                    else:
                        # 0으로 2컬럼 만들어서 채워줘야함.
                        l_NULL = ['0' for i in range(0, len(ped_Sample_id))]
                        sr_NULL = pd.concat([pd.Series(l_NULL), pd.Series(l_NULL)], axis=1)
                        sr_NULL.index = pd.Index(ped_Sample_id)

                        l_HIBAG.append(sr_NULL)

                df_HIBAG = pd.concat(l_HIBAG, axis=1)


                ped_Family_ID = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Paternal_ID = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Maternal_ID = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Sex = ['0' for i in range(0, len(ped_Sample_id))]
                ped_Affection = ['0' for i in range(0, len(ped_Sample_id))]


                df_HIBAG.index = pd.MultiIndex.from_arrays([ped_Family_ID, ped_Sample_id, ped_Paternal_ID, ped_Maternal_ID, ped_Sex, ped_Affection])

                print(df_HIBAG.head())

                df_HIBAG.to_csv(_OUTPUT+'.hibag.ped', sep='\t', header=False, index=True)

            else:
                # 찾은 HIBAG 파일이 하나도 없는 경우
                print("No HIBAG file for any HLA! Please check input file directory and files again.")
                sys.exit()



        else:
            print("\n\nGiven _PLATFORM argument is worng. Check again!\n")
            sys.exit()






if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################
    #
    # Created by Wansun Choi, 2017/12/06
    # CONVERTER.py
    #
    # SPECIFICATION : This program converts results files from HLA Imputation software 'HIBAG' or 'Axiom HLA Typing Analysis Services' to .ped file that can be used in SNP2HLA
    #
    # INPUT file : Result Files from 'HIBAG(https://bioconductor.org/packages/release/bioc/html/HIBAG.html)' or 'Axiom HLA Typing Analysis Services(https://www.thermofisher.com/order/catalog/product/000.911)'
    # OUTPUT file : .ped file
    #
    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
    parser.add_argument("-p", "--PLATFORM", help="\nChoose a platform [HIBAG] or [AXIOM]\n(HIBAG - https://bioconductor.org/packages/release/bioc/html/HIBAG.html)\n(Axiom HLA Typing Analysis Services - https://www.thermofisher.com/order/catalog/product/000.911)\n\n",
                        choices=['AXIOM', 'HIBAG'], nargs='?', metavar="platform", required=True)
    parser.add_argument("-i", help="\nyour input file(s) prefix\n\n", nargs='+', metavar='filename', required=True)
    parser.add_argument("-o", "--OUTPUT", help="\nyour output file prefix\n\n", nargs='?', metavar='prefix', default='CONVERTED')



    ##### <for Test> #####

    # args = parser.parse_args(["-p", "AXIOM", "-i", "./data_CONVERTER/AxiomHLA_TestResult", "-o", "TEST_hg19.converted.AXIOM"]) # 폴더

    # args = parser.parse_args(["-p", "HIBAG", "-i", "./data_CONVERTER/HIBAG_TestResult_HLA-A.txt", "-o", "TEST_hg19.converted.HIBAG"]) # single file
    # args = parser.parse_args(["-p", "HIBAG", "-i", "./data_CONVERTER/HIBAG_TestResult_HLA-A.txt", "HIBAG_TestResult_HLA-B.txt", "-o", "TEST_hg19.converted.HIBAG"]) # multiple files
    # args = parser.parse_args(["-p", "HIBAG", "-i", "./data_CONVERTER/HIBAG_TestResult_HLA-A.txt", "HIBAG_TestResult_HLA-B.txt", "./asdfsd/fadsfas/HIBAG_TestResult_HLA-C.txt", "./dfadsf/adsfdasf/dsfadfd/HIBAG_TestResult_HLA-DPB1.txt", "-o", "TEST_hg19.converted.HIBAG"]) # multiple files




    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)

    # print(parser.print_help())

    print("Parameters are given.")
    # argument passing

    # main function execution
    CONVERTER(args.PLATFORM, args.i, args.OUTPUT)

