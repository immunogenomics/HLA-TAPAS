# -*- coding: utf-8 -*-

import os, sys, re
import subprocess

from MakeReference.src.HLAtoSequences import HLAtoSequences
from MakeReference.src.encodeVariants import encodeVariants
from MakeReference.src.encodeHLA import encodeHLA

from MakeReference.src.ATtrick import ATtrick
from MakeReference.src.redefineBPv1BH import redefineBP

########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def MakeReference_v2(_CHPED, _OUT, _hg, _dictionary_AA, _dictionary_SNPS, _variants=None,
                     _java_heap_mem='2g', _java_stack_mem='1g',
                     _p_dependency="dependency/", f_save_intermediates=False, f_phasing=False, _nthreads=1):


    ########## < Core Variables > ##########

    ### [1] Major Path Variables
    if os.path.exists(_p_dependency):
        p_dependency = _p_dependency
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Folder for dependent software('{}') can't be found. Please check '--dependency' argument again.".format(_p_dependency))
        sys.exit()

    # dependent software.
    _p_plink = os.path.join(p_dependency, "plink")
    _p_beagle = os.path.join(p_dependency, "beagle.jar")
    _p_linkage2beagle = os.path.join(p_dependency, "linkage2beagle.jar")
    _p_beagle2vcf = os.path.join(p_dependency, "beagle2vcf.jar")

    if not os.path.exists(_p_plink):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(_p_dependency))
        sys.exit()

    if not os.path.exists(_p_beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'Beagle 4.1' (https://faculty.washington.edu/browning/beagle/b4_1.html#download) in '{0}'\n".format(_p_dependency))
        sys.exit()

    if not os.path.exists(_p_linkage2beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'linkage2beagle.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) in '{0}'\n".format(_p_dependency))
        sys.exit()

    if not os.path.exists(_p_beagle2vcf):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'beagle2vcf.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) in '{0}'\n".format(_p_dependency))
        sys.exit()



    ### [2] Memory representation check.
    p_Mb = re.compile(r'\d+m')
    p_Gb = re.compile(r'\d+[gG]')

    if not (bool(p_Mb.match(_java_heap_mem)) or bool(p_Gb.match(_java_heap_mem))):
        print(std_ERROR_MAIN_PROCESS_NAME + "Given Java memory value('{}') has bizzare representation. Please check '--mem' argument again.".format(_java_heap_mem))
        sys.exit()



    ### [3] Dictionary Files

    _dictionary_AA_seq = _dictionary_AA + ".txt" # From now on, official extension of HLA sequence information dictionary is ".txt". (2018. 9. 25.)
    _dictionary_AA_map = _dictionary_AA + ".map"

    _dictionary_SNPS_seq = _dictionary_SNPS + ".txt"
    _dictionary_SNPS_map = _dictionary_SNPS + ".map"


    if not os.path.exists(_dictionary_AA_map):
        print(std_ERROR_MAIN_PROCESS_NAME + "AA dictionary map file can't be found('{0}'). Please check '--dict-AA' argument again.\n".format(_dictionary_AA_map))
        sys.exit()

    if not os.path.exists(_dictionary_AA_seq):
        print(std_ERROR_MAIN_PROCESS_NAME + "AA dictionary seq file can't be found('{0}'). Please check '--dict-AA' argument again.\n".format(_dictionary_AA_seq))
        sys.exit()

    if not os.path.exists(_dictionary_SNPS_map):
        print(std_ERROR_MAIN_PROCESS_NAME + "SNPS dictionary map file can't be found('{0}'). Please check '--dict-SNPS' argument again.\n".format(_dictionary_SNPS_map))
        sys.exit()

    if not os.path.exists(_dictionary_SNPS_seq):
        print(std_ERROR_MAIN_PROCESS_NAME + "SNPS dictionary seq file can't be found('{0}'). Please check '--dict-SNPS' argument again.\n".format(_dictionary_SNPS_seq))
        sys.exit()



    ### Intermediate path.
    OUTPUT = _OUT if not _OUT.endswith('/') else _OUT.rstrip('/')
    if bool(os.path.dirname(OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)



    ########## < Core Variables 2 > ##########

    ### Flag for plain SNP markers.
    f_plain_SNP = bool(_variants)

    # Input 1 : HLA type data
    HLA_DATA = _CHPED

    if not os.path.exists(_CHPED):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA type can't be found('{}'). Please check '--chped' argument again.".format(_CHPED))
        sys.exit()

    # Input 2 : Plain SNP data
    if f_plain_SNP:

        if not os.path.exists(_variants+'.bed'):
            print(std_ERROR_MAIN_PROCESS_NAME + "One of SNP data can't be found('{}'). Please check '--variants' argument again.".format(_variants+'.bed'))
            sys.exit()
        if not os.path.exists(_variants+'.bim'):
            print(std_ERROR_MAIN_PROCESS_NAME + "One of SNP data can't be found('{}'). Please check '--variants' argument again.".format(_variants+'.bim'))
            sys.exit()
        if not os.path.exists(_variants+'.fam'):
            print(std_ERROR_MAIN_PROCESS_NAME + "One of SNP data can't be found('{}'). Please check '--variants' argument again.".format(_variants+'.fam'))
            sys.exit()

        SNP_DATA = _variants
        SNP_DATA2 = os.path.join(INTERMEDIATE_PATH, os.path.basename(_variants))

    else:
        print(std_WARNING_MAIN_PROCESS_NAME + "SNP data wasn't given. MakeReference_v2 will generate a reference panel only with HLA type information.")


    # Output prefix
    AA_CODED = OUTPUT + '.AA.CODED'
    HLA_CODED = OUTPUT + ".HLA"
    SNPS_CODED = OUTPUT + '.SNPS.CODED'

    plink = ' '.join([_p_plink, "--noweb", "--silent"])
    beagle = ' '.join(["java", "-Xmx{}".format(_java_heap_mem), "-Xss{}".format(_java_stack_mem), "-jar", _p_beagle])
    linkage2beagle = ' '.join(["java", "-Xmx{}".format(_java_heap_mem), "-jar", _p_linkage2beagle])
    beagle2vcf = ' '.join(["java", "-Xmx{}".format(_java_heap_mem), "-jar", _p_beagle2vcf])




    ########## <Flags for Code Block> ##########

    ENCODE_AA = 1
    ENCODE_HLA = 1
    ENCODE_SNPS = 1

    EXTRACT_FOUNDERS = 1
    MERGE = 1
    QC = 1

    PREPARE = 1
    PHASE = 1
    CLEANUP = 0

    # (2019. 01. 10.) Last three code blocks won't be implemented



    ########## <Making Reference Panel> ##########

    print(std_MAIN_PROCESS_NAME + "Making Reference Panel for \"{0}\"".format(OUTPUT))
    index = 1

    if ENCODE_AA:

        '''
        echo "[$i] Generating amino acid sequences from HLA types.";  @ i++
        ./HLAtoSequences.pl $HLA_DATA HLA_DICTIONARY_AA.txt AA > $OUTPUT.AA.ped
        cp HLA_DICTIONARY_AA_hg19.map $OUTPUT.AA.map # hg19
        # cp HLA_DICTIONARY_AA.map $OUTPUT.AA.map

        echo "[$i] Encoding amino acids positions." ;  @ i++
        ./encodeVariants.pl $OUTPUT.AA.ped $OUTPUT.AA.map $OUTPUT.AA.CODED

        plink --file $OUTPUT.AA.CODED --missing-genotype 0 --make-bed --out $OUTPUT.AA.TMP
        awk '{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}' $OUTPUT.AA.TMP.bim | grep -v INS | cut -f2 > to_remove
        plink --bfile $OUTPUT.AA.TMP --exclude to_remove --make-bed --out $OUTPUT.AA.CODED

        # rm $OUTPUT.AA.TMP*; rm to_remove
        # rm $OUTPUT.AA.???
        '''

        print("[{}] Generating Amino acid(AA)sequences from HLA types.".format(index))

        ### (1) Sequence File ( *.AA.{ped,map} ) ###
        HLAtoSequences(HLA_DATA, _dictionary_AA_seq, "AA", _out=OUTPUT)
        os.system(' '.join(["cp", _dictionary_AA_map, OUTPUT + '.AA.map']))

        index += 1


        print("[{}] Encoding Amino acids positions.".format(index))

        ### (2) pre-CODED File ( *.AA.CODED.{ped,map},  *.AA.CODED.factors ) ###
        encodeVariants(OUTPUT + '.AA.ped', OUTPUT + '.AA.map', OUTPUT + '.AA.CODED')  # previously "enCODED".

        index += 1


        ### (3) stuffs to remove ( *.AA.TMP.{bed,bim,fam}, to_remove ) ###
        command = ' '.join(
            [plink, "--file", OUTPUT + '.AA.CODED', "--missing-genotype 0", "--make-bed", "--out", OUTPUT + '.AA.TMP'])
        # print(command)
        os.system(command)

        command = ' '.join(
            ["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.AA.TMP.bim', "|",
             "cut -f2", ">",
             os.path.join(INTERMEDIATE_PATH, "to_remove")])

        """
        In the previous framework which was created by Sherman Jia, The insertions were dealt as a marker "INS".
        In the new version of Framework, marker label is "INDEL".
        """

        # print(command)
        os.system(command)


        ### (4) Final Encoded outputs ( *.AA.CODED.{bed,bim,fam,nosex,log} ) ###

        command = ' '.join(
            [plink, "--bfile", OUTPUT + '.AA.TMP',
             "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"),
             "--make-bed",
             "--out", OUTPUT + '.AA.CODED'])
        # print(command)
        os.system(command)


        """
        Generated outputs :
            (1) Sequence File ( *.AA.{ped,map} )
            (2) pre-CODED File ( *.AA.CODED.{ped,map}, *.AA.CODED.factors )
            (3) stuffs to remove ( *.AA.TMP.{bed,bim,fam}, to_remove )
            (4) Final Encoded outputs ( *.AA.CODED.{bed,bim,fam,nosex,log}
            
            
        Final outputs :
            - *.AA.CODED.{bed,bim,fam,factors,nosex,log}
        
        Outputs to remove :
            - *.AA.{ped,map}
            - *.AA.TMP.*
            - to_remove
            - *.AA.CODED.{ped,map}
            
        """

        if not f_save_intermediates:

            os.system("rm " + (OUTPUT + ".AA.ped"))
            os.system("rm " + (OUTPUT + ".AA.map"))
            os.system("rm " + (OUTPUT + ".AA.TMP.*"))
            os.system("rm " + os.path.join(INTERMEDIATE_PATH, "to_remove"))
            os.system("rm " + OUTPUT + ".AA.CODED.ped")
            os.system("rm " + OUTPUT + ".AA.CODED.map")
            os.system("rm " + OUTPUT + ".AA.CODED.factors")



    if ENCODE_HLA:

        print("[{}] Encoding HLA alleles.".format(index))

        ### (1) Encoded HLA ( *.HLA.{ped,map} ) ###
        encodeHLA(HLA_DATA, OUTPUT + ".HLA", _hg)

        ### (2) Final Encoded Outputs ( *.HLA.{bed,bim,fam,nosex,log} ) ###
        command = ' '.join([plink, "--file", OUTPUT + '.HLA', "--make-bed", "--out", OUTPUT + '.HLA'])
        # print(command)
        os.system(command)

        index += 1


        """
        Generaed outputs :
            (1) Encoded HLA ( *.HLA.{ped,map} )
            (2) Final Encoded Outputs ( *.HLA.{bed,bim,fam}, *.HLA.{nosex,log} )
            
        Final outputs :
            - *.HLA.{bed,bim,fam,factors,nosex,log}
            
        Outputs to remove :
            - *.HLA.{ped,map}
        """

        if not f_save_intermediates:
            os.system("rm " + (OUTPUT + ".HLA.ped"))
            os.system("rm " + (OUTPUT + ".HLA.map"))




    if ENCODE_SNPS:

        print("[{}] Generating DNA(SNPS) sequences from HLA types.".format(index))

        ### (1) Sequence File ( *.SNPS.{ped,map} )
        HLAtoSequences(HLA_DATA, _dictionary_SNPS_seq, "SNPS", OUTPUT)

        command = ' '.join(["cp", _dictionary_SNPS_map, OUTPUT + '.SNPS.map'])
        # print(command)
        os.system(command)

        index += 1


        print("[{}] Encoding SNP positions.".format(index))

        ### (2) pre-CODED File ( *.SNPS.CODED.{ped,map} )
        encodeVariants(OUTPUT + '.SNPS.ped', OUTPUT + '.SNPS.map', OUTPUT + '.SNPS.CODED')

        index += 1


        ### (3) Stuffs to remove ( *.SNPS.TMP.{bed,bim,fam}, to_remove )
        command = ' '.join([plink, "--file", OUTPUT + '.SNPS.CODED', "--missing-genotype 0", "--make-bed", "--out",
                            OUTPUT + '.SNPS.TMP'])
        # print(command)
        os.system(command)


        command = ' '.join(
            ["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.SNPS.TMP.bim', "|", "cut -f2", ">",
             os.path.join(INTERMEDIATE_PATH, "to_remove")])
        # print(command)
        os.system(command)


        ### (4) Final Encoded outputs ( *.SNPS.CODED.{bed,bim,fam,nosex,log} )
        command = ' '.join(
            [plink, "--bfile", OUTPUT + '.SNPS.TMP',
             "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"),
             "--make-bed",
             "--out", OUTPUT + '.SNPS.CODED'])
        # print(command)
        os.system(command)


        """
        Generated outputs :
            (1) Sequence File ( *.SNPS.{ped,map} )
            (2) pre-CODED File ( *.SNPS.CODED.{ped,map}, *.SNPS.CODED.factors )
            (3) stuffs to remove ( *.SNPS.TMP.{bed,bim,fam}, to_remove )
            (4) Final Encoded outputs ( *.SNPS.CODED.{bed,bim,fam,nosex,log}
            
                      
        Final outputs :
            - *.SNPS.CODED.{bed,bim,fam,factors,nosex,log}
        
        Outputs to remove :
            - *.SNPS.{ped,map}
            - *.SNPS.TMP.*
            - to_remove
            - *.SNPS.CODED.{ped,map}
                    
        """

        if not f_save_intermediates:

            os.system("rm " + (OUTPUT + ".SNPS.ped"))
            os.system("rm " + (OUTPUT + ".SNPS.map"))
            os.system("rm " + (OUTPUT + ".SNPS.TMP.*"))
            os.system("rm " + os.path.join(INTERMEDIATE_PATH, "to_remove"))
            os.system("rm " + (OUTPUT + ".SNPS.CODED.ped"))
            os.system("rm " + (OUTPUT + ".SNPS.CODED.map"))
            os.system("rm " + (OUTPUT + ".SNPS.CODED.factors"))




    # Above 3 code blocks are implmented no matter `_variants`(Normal genomewide SNPs) is given or not.
    # Next parts are implented sligtly differently depending on `_variants`.


    if _variants:

        ### `_variants` is given.
        # Almost same as original "MakeReference".

        if EXTRACT_FOUNDERS:

            print("[{}] Extracting founders.".format(index))

            """
            if ($EXTRACT_FOUNDERS) then
                echo "[$i] Extracting founders."; @ i++
                plink --bfile $SNP_DATA --filter-founders --mind 0.3 --alleleACGT --make-bed --out $SNP_DATA.FOUNDERS
            
                # Initial QC on Reference SNP panel
                plink --bfile $SNP_DATA.FOUNDERS --hardy        --out $SNP_DATA.FOUNDERS.hardy  # 진짜 92명에 대해 position별로 HWE test한 결과
                plink --bfile $SNP_DATA.FOUNDERS --freq         --out $SNP_DATA.FOUNDERS.freq   # 실제 --freq 옵션이 allele frequency계산해주는 옵션임. 
                plink --bfile $SNP_DATA.FOUNDERS --missing      --out $SNP_DATA.FOUNDERS.missing
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.hardy.hwe      | awk ' $9 < 0.000001 { print $2 }' | sort -u > remove.snps.hardy 
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.freq.frq       | awk ' $5 < 0.01 { print $2 } '             > remove.snps.freq
                awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.missing.lmiss  | awk ' $5 > 0.05 { print $2 } '              > remove.snps.missing
                cat remove.snps.*                                            | sort -u                                     > all.remove.snps
            
                plink --bfile $SNP_DATA.FOUNDERS --allow-no-sex --exclude all.remove.snps --make-bed --out $SNP_DATA.FOUNDERS.QC
            
                # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file
            
                plink --bfile $OUTPUT.HLA --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.HLA.FOUNDERS
                plink --bfile $OUTPUT.SNPS.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.SNPS.FOUNDERS
                plink --bfile $OUTPUT.AA.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.AA.FOUNDERS
            
                rm remove.snps.*
            endif
            """

            ### (1) --filter-founders to variants data ( *.FOUNDERS )
            command = ' '.join([plink, "--bfile", SNP_DATA, "--filter-founders", "--mind 0.3", "--alleleACGT", "--make-bed", "--out", SNP_DATA2+'.FOUNDERS'])
            # print(command)
            os.system(command)


            ### (2) QC ( *.FOUNDERS.{hardy,freq,missing )
            # Initial QC on Reference SNP panel
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--hardy", "--out", SNP_DATA2+'.FOUNDERS.hardy'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--freq", "--out", SNP_DATA2+'.FOUNDERS.freq'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--missing", "--out", SNP_DATA2+'.FOUNDERS.missing'])
            # print(command)
            os.system(command)


            ### (3) Stuffs to remove ( remove.snps.hardy, remove.snps.freq, remove.snps.missing, all.remove.snps )
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.hardy.hwe', "|", "awk", "' $9 < 0.000001 { print $2 }'", "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.hardy")])
            # print(command)
            os.system(command)
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.freq.frq', "|", "awk", "' $5 < 0.01 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.freq")])
            # print(command)
            os.system(command)
            command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.missing.lmiss', "|", "awk", "' $5 > 0.05 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.missing")])
            # print(command)
            os.system(command)
            command = ' '.join(["cat", os.path.join(INTERMEDIATE_PATH, "remove.snps.*"), "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "all.remove.snps")])
            # print(command)
            os.system(command)


            ### (4) Filtering out Quality-controled FOUNDERS ( *.FOUNDERS.QC )
            command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--allow-no-sex", "--exclude", os.path.join(INTERMEDIATE_PATH, "all.remove.snps"), "--make-bed", "--out", SNP_DATA2+'.FOUNDERS.QC'])
            # print(command)
            os.system(command)

            # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

            ### (5) --filter-founders to HLA information ( *.{HLA,AA,SNPS}.FOUNDERS )
            command = ' '.join([plink, "--bfile", OUTPUT+'.HLA', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.HLA.FOUNDERS'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", OUTPUT+'.SNPS.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.SNPS.FOUNDERS'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", OUTPUT+'.AA.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.AA.FOUNDERS'])
            # print(command)
            os.system(command)


            """
            Generated outputs :
                (1) --filter-founders to variants data ( *.FOUNDERS )
                (2) QC ( *.FOUNDERS.{hardy,freq,missing )
                (3) Stuffs to remove ( remove.snps.hardy, remove.snps.freq, remove.snps.missing, all.remove.snps )
                (4) Filtering out Quality-controled FOUNDERS ( *.FOUNDERS.QC )
                (5) --filter-founders to HLA information ( *.{HLA,AA,SNPS}.FOUNDERS )
                
            Final outputs :
                - *.FOUNDERS.QC
                - *.{HLA,AA,SNPS}.FOUNDERS
                
            Outputs to remove :
                - *.FOUNDERS
                - *.FOUNDERS.{hardy,freq,missing}
                - remove.snps.hardy, remove.snps.freq, remove.snps.missing, all.remove.snps
                
            """


            if not f_save_intermediates:

                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.bed"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.bim"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.fam"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.log"))
                if os.path.exists(SNP_DATA2 + ".FOUNDERS.nosex"):
                    os.system("rm " + (SNP_DATA2 + ".FOUNDERS.nosex"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.hardy.*"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.freq.*"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.missing.*"))
                os.system("rm " + os.path.join(INTERMEDIATE_PATH, "remove.snps.*"))
                os.system("rm " + os.path.join(INTERMEDIATE_PATH, "all.remove.snps"))


            index += 1


        if MERGE:

            print("[{}] Merging SNP, HLA, and amino acid datasets.".format(index))

            """
            echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
            echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
            echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
            echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
            plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
            rm $OUTPUT.HLA.???
            rm $OUTPUT.AA.CODED.???
            rm $OUTPUT.SNPS.CODED.???
            rm merge_list
    
            """

            ### (1) Stuffs to merge ( merge_list )
            TMP_merged_list = os.path.join(INTERMEDIATE_PATH, "merge_list")


            command = ' '.join(["echo", OUTPUT + '.HLA.FOUNDERS.bed', OUTPUT + '.HLA.FOUNDERS.bim', OUTPUT + '.HLA.FOUNDERS.fam', ">", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(["echo", OUTPUT + '.AA.FOUNDERS.bed', OUTPUT + '.AA.FOUNDERS.bim', OUTPUT + '.AA.FOUNDERS.fam', ">>", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(["echo", OUTPUT + '.SNPS.FOUNDERS.bed', OUTPUT + '.SNPS.FOUNDERS.bim', OUTPUT + '.SNPS.FOUNDERS.fam', ">>", TMP_merged_list])
            # print(command)
            os.system(command)


            ### (2) Merging the above stuffs ( *.MERGED.FOUNDERS )
            command = ' '.join(
                [plink, "--bfile", SNP_DATA2 + '.FOUNDERS.QC', "--merge-list", TMP_merged_list, "--make-bed", "--out",
                 OUTPUT + '.MERGED.FOUNDERS'])
            # print(command)
            os.system(command)


            """
            Generated Outputs :
                (1) Stuffs to merge ( merge_list )
                (2) Merging the above stuffs ( *.MERGED.FOUNDERS )
            
            Final outputs : 
                - *.MERGED.FOUNDERS
            
            Outputs to remove :
                - *.{AA,HLA,SNPS}.FOUNDERS.{bed,bim,fam}
                - merge_list
            
            
            """

            if not f_save_intermediates:

                os.system("rm " + (OUTPUT + ".HLA.*"))
                os.system("rm " + (OUTPUT + ".AA.*"))
                os.system("rm " + (OUTPUT + ".SNPS.*"))
                os.system("rm " + (SNP_DATA2 + ".FOUNDERS.QC.*"))
                os.system("rm " + TMP_merged_list)



            index += 1


        if QC:

            print("[{}] Performing quality control.".format(index))

            """
            plink --bfile $OUTPUT.MERGED.FOUNDERS --freq --out $OUTPUT.MERGED.FOUNDERS.FRQ
            awk '{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > all.remove.snps
            awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order
    
            # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
            plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --exclude all.remove.snps --geno 0.5 --make-bed --out $OUTPUT
    
            # Calculate allele frequencies
            plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
            rm $SNP_DATA.FOUNDERS.*
            rm $OUTPUT.MERGED.FOUNDERS.*
            rm $OUTPUT.*.FOUNDERS.???
            rm allele.order
            rm all.remove.snps
    
            """

            TMP_all_remove_snps = os.path.join(INTERMEDIATE_PATH, "all.remove.snps")
            # TMP_allele_order = os.path.join(INTERMEDIATE_PATH, "allele.order")
            TMP_allele_order = OUTPUT + ".refallele"


            ### (1) Frequency file to use in filtering out some snps ( *.MERGED.FOUNDERS.FRQ )
            command = ' '.join(
                [plink, "--bfile", OUTPUT + '.MERGED.FOUNDERS', "--freq", "--out", OUTPUT + '.MERGED.FOUNDERS.FRQ'])
            # print(command)
            os.system(command)

            ### (2) List up SNPs which have extreme allele frequency. ( all.remove.snps )
            command = ' '.join(["awk", "'{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}'",
                                OUTPUT + '.MERGED.FOUNDERS.FRQ.frq', ">", TMP_all_remove_snps])
            # print(command)
            os.system(command)

            ### (3) Making reference allele ( allele.order )
            # command = ' '.join(["awk", '\'{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}\'', OUTPUT+'.MERGED.FOUNDERS.FRQ.frq', ">", TMP_allele_order])
            command = ' '.join(
                ["awk", '\'{if (NR > 1){if (($3 == "a" && $4 == "p") || ($4 == "a" && $3 == "p")){print $2 "\tp"}}}\'',
                 OUTPUT + '.MERGED.FOUNDERS.FRQ.frq', ">", TMP_allele_order])
            # print(command)
            os.system(command)


            ### (4) Filtering out SNPs which have extreme allele frequency(The reference panel for association test.) ( *.{bed,bim,fam} )
            # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
            command = ' '.join(
                [plink, "--bfile", OUTPUT + '.MERGED.FOUNDERS', "--a1-allele", TMP_allele_order, "--exclude",
                 TMP_all_remove_snps, "--geno 0.5", "--make-bed", "--out", OUTPUT])
            # print(command)
            os.system(command)


            ### (5) Allele frequency info. of output reference panel ( *.FRQ )
            # Calculate allele frequencies
            command = ' '.join([plink, "--bfile", OUTPUT, "--keep-allele-order", "--freq", "--out", OUTPUT + '.FRQ',
                                "--a1-allele", TMP_allele_order])
            # print(command)
            os.system(command)




            """
            Generated Outputs :
                (1) Frequency file to use in filtering out some snps ( *.MERGED.FOUNDERS.FRQ )
                (2) List up SNPs which have extreme allele frequency. ( all.remove.snps )
                (3) Making reference allele ( allele.order )
                (4) Filtering out SNPs which have extreme allele frequency(The reference panel for association test.) ( *.{bed,bim,fam} )
                (5) Allele frequency info. of output reference panel ( *.FRQ )

            Final outputs : 
                - *.{bed,bim,fam}
                - *.FRQ

            Outputs to remove :
                - *.MERGED.FOUNDERS.FRQ
                - all.remove.snps
                - allele.order
                
            """

            if not f_save_intermediates:

                os.system("rm " + (OUTPUT + ".MERGED.FOUNDERS.*"))
                os.system("rm " + (OUTPUT + ".FRQ.log"))
                os.system("rm " + (OUTPUT + ".log"))
                os.system("rm " + TMP_all_remove_snps)
                os.system("rm " + TMP_allele_order)

                if os.path.exists(OUTPUT + ".FRQ.nosex"):
                    os.system("rm " + (OUTPUT + ".FRQ.nosex"))
                if os.path.exists(OUTPUT + ".nosex"):
                    os.system("rm " + (OUTPUT + ".nosex"))

            index += 1



        if f_phasing:
    

            if PREPARE:

                print("[{}] Preparing files for Beagle.".format(index))

                """
                [Source from Buhm Han.]
    
                awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
                plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
                awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
                cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
    
                echo "[$i] Converting to beagle format.";  @ i++
                linkage2beagle pedigree=$OUTPUT.nopheno.ped data=$OUTPUT.dat beagle=$OUTPUT.bgl standard=true > $OUTPUT.bgl.log
    
    
                [Source from Yang.]
    
                awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
                plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
                plink --bfile $OUTPUT --recode --transpose --out $OUTPUT
                # awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
                # cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
    
                echo "[$i] Converting to beagle format.";  @ i++
                beagle2vcf -fnmarker $OUTPUT.markers -fnfam $OUTPUT.fam -fngtype $OUTPUT.tped -fnout $OUTPUT.vcf
    
                I will make this code block based on source given by Yang. for now.
    
                """

                # ATtrick : 'P' to 'T', 'A' to 'A'
                [bim_ATtrick, a1_allele_ATtrick] = ATtrick(OUTPUT + '.bim', OUTPUT)


                # command = ' '.join(["awk", '\'{print $2" "$4" "$5" "$6}\'', OUTPUT + '.bim', ">", OUTPUT + '.markers'])
                command = ' '.join(["awk", '\'{print $2" "$4" "$6" "$5}\'', bim_ATtrick, ">", OUTPUT + '.ATtrick.markers'])
                # print(command)
                os.system(command)
                """
                Plink works by setting ALT allele as a1-allele, which is the 5th column of bim file.
                However, VCF file sets a2-allele, which is the 6th column of plink bim file, as ALT allele.
                
                beagle2vcf converts 4th column of *.markers file to ALT allele of newly generated vcf file.
                That's why I reordered ($5, $6) to ($6, $5).
                """

                # Manipulate duplicated Base poistion
                redefined_markers = redefineBP(OUTPUT + '.ATtrick.markers', OUTPUT + '.markers')

                # Applying the above manipulated base position information.
                # command = [
                #     "paste <(cut -f1-3 %s) <(awk '{print $2}' %s) <(cut -f5- %s) > %s" % (bim_ATtrick, redefined_markers, bim_ATtrick, OUTPUT+'.ATtrick.redefined.bim')
                # ]
                # print(command)
                # # os.system(command)
                # subprocess.call(command)

                command = ' '.join(
                    [plink, "--bed", OUTPUT+'.bed', "--bim", bim_ATtrick, "--fam", OUTPUT+'.fam',
                     "--keep-allele-order", "--recode", "--alleleACGT", "--out", OUTPUT+'.ATtrick',
                     "--a1-allele", a1_allele_ATtrick])
                # print(command)
                os.system(command)

                command = ' '.join(["awk", '\'{print "M " $2}\'', OUTPUT + '.ATtrick.map', ">", OUTPUT + '.ATtrick.dat'])
                # print(command)
                os.system(command)

                command = ' '.join(["cut -d ' ' -f1-5,7-", OUTPUT + '.ATtrick.ped', ">", OUTPUT + '.ATtrick.nopheno.ped'])
                # print(command)
                os.system(command)

                index += 1

                print("[{}] Converting PLINK to BEAGLE format.".format(index))

                command = ' '.join([linkage2beagle, "pedigree=" + OUTPUT + '.ATtrick.nopheno.ped', "data=" + OUTPUT + '.ATtrick.dat',
                                    "beagle=" + OUTPUT + '.ATtrick.bgl', "standard=true", ">", OUTPUT + '.ATtrick.bgl.log'])
                # print(command)
                os.system(command)

                index += 1

                # for Beagle 4.1
                print("[{}] Converting BEAGLE to VCF format.".format(index))

                command = ' '.join([beagle2vcf, '6', redefined_markers, OUTPUT + '.ATtrick.bgl', '0', '>', OUTPUT+'.bgl.vcf'])
                # print(command)
                os.system(command)

                index += 1


                if not f_save_intermediates:
                    # os.system("rm {}".format())
                    os.system("rm {}".format(OUTPUT + '.ATtrick.bgl.log'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.bgl'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.log'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.dat'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.nopheno.ped'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.map'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.ped'))
                    os.system("rm {}".format(OUTPUT + '.ATtrick.markers'))
                    os.system("rm {}".format(bim_ATtrick))
                    os.system("rm {}".format(a1_allele_ATtrick))



            if PHASE:

                print("[{}] Phasing reference using Beagle4.1.".format(index))

                '''
                # Beagle v3.x.x
                beagle unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log
    
                # Beagle v4.1
                beagle gt=OUTPUT+'.ATtrick.redefined.vcf' out=OUTPUT+'.ATtrick.redefined.phased' nthreads=1 impute=false niterations=10 lowmem=true >> $OUTPUT.bgl.log
    
                '''

                # command = ' '.join([beagle, "gt={}".format(OUTPUT+'.bglv4.bgl.vcf'),
                #                     "nthreads=1", "impute=false",
                #                     "niterations=10", "lowmem=true", "out={}".format(OUTPUT+'.bglv4.bgl.phased'),
                #                     '>', OUTPUT+'.bglv4.bgl.phased.vcf.log'])

                command = ' '.join([beagle, "gt={}".format(OUTPUT+'.bgl.vcf'),
                                    "impute=false", "nthreads={}".format(_nthreads),
                                    "niterations=5", "lowmem=true", "out={}".format(OUTPUT+'.bgl.phased')])
                # print(command)

                try:
                    # os.system(command)
                    f_log = open(OUTPUT+'.bgl.phased.vcf.log', 'w')
                    subprocess.run(re.split(r'\s+',command), check=True, stdout=f_log, stderr=f_log)

                except subprocess.CalledProcessError:
                    # fail.
                    print(std_ERROR_MAIN_PROCESS_NAME + "Phasing failed. See log file('{}').".format(OUTPUT+'.bgl.phased.vcf.log'))
                    sys.exit()
                else:
                    # succeed.
                    f_log.close()
                    if not f_save_intermediates:
                        # os.system("rm {}".format())
                        os.system("rm {}".format(OUTPUT + '.bgl.vcf'))

                        # remove redundant log file.
                        os.system("rm {}".format(OUTPUT+'.bgl.phased.log'))

                index += 1



        if CLEANUP:

            print("[{}] Removing unnecessary files.".format(index))

            '''
            rm $OUTPUT.nopheno.ped
            rm $OUTPUT.bgl.gprobs
            rm $OUTPUT.bgl.r2
            rm $OUTPUT.bgl
            rm $OUTPUT.ped
            rm $OUTPUT.map
            rm $OUTPUT.dat
            rm $OUTPUT.phasing.log
            '''

            rm_tlist = ('.nopheno.ped', '.bgl.gprobs', '.bgl.r2', '.bgl', '.ped', '.map', '.dat')

            for i in rm_tlist:
                print("rm " + OUTPUT + i)
                os.system("rm " + OUTPUT + i)

            index += 1


        print("[{}] Making reference panel for HLA-AA,SNPS,HLA and Normal variants(SNPs) is Done!".format(index))

        __return__ = OUTPUT



    else:

        ### No `_variants` is given.
        # Only for HLA markers.
        # "EXTRACT_FOUNDERS" code block is not included.


        __HLA__ = None
        __AA__ = None
        __SNPS__ = None
        __MERGED__ = None


        if EXTRACT_FOUNDERS:

            print("[{}] Extracting founders.".format(index))


            __HLA__ = OUTPUT + ".HLA"
            __AA__ = OUTPUT + ".AA.CODED"
            __SNPS__ = OUTPUT + ".SNPS.CODED"


            command = ' '.join([plink, "--bfile", __HLA__, "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.HLA.FOUNDERS'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", __SNPS__, "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.SNPS.FOUNDERS'])
            # print(command)
            os.system(command)
            command = ' '.join([plink, "--bfile", __AA__, "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.AA.FOUNDERS'])
            # print(command)
            os.system(command)


            if not f_save_intermediates:

                os.system("rm " + OUTPUT+".HLA.bed")
                os.system("rm " + OUTPUT+".HLA.bim")
                os.system("rm " + OUTPUT+".HLA.fam")
                os.system("rm " + OUTPUT+".HLA.log")

                os.system("rm " + OUTPUT+".SNPS.CODED.bed")
                os.system("rm " + OUTPUT+".SNPS.CODED.bim")
                os.system("rm " + OUTPUT+".SNPS.CODED.fam")
                os.system("rm " + OUTPUT+".SNPS.CODED.log")

                os.system("rm " + OUTPUT+".AA.CODED.bed")
                os.system("rm " + OUTPUT+".AA.CODED.bim")
                os.system("rm " + OUTPUT+".AA.CODED.fam")
                os.system("rm " + OUTPUT+".AA.CODED.log")

                if os.path.exists(__HLA__+".nosex"):
                    os.system("rm " + __HLA__ + ".nosex")
                if os.path.exists(__SNPS__+".nosex"):
                    os.system("rm " + __SNPS__ + ".nosex")
                if os.path.exists(__AA__+".nosex"):
                    os.system("rm " + __AA__ + ".nosex")



        if MERGE:

            print("[{}] Merging SNP, HLA, and amino acid datasets.".format(index))

            """
            echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
            echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
            echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
            echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
            plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
            rm $OUTPUT.HLA.???
            rm $OUTPUT.AA.CODED.???
            rm $OUTPUT.SNPS.CODED.???
            rm merge_list
    
            """

            # Output of previous code block.
            __AA__ = OUTPUT+'.AA.FOUNDERS'
            __HLA__ = OUTPUT+'.HLA.FOUNDERS'
            __SNPS__ = OUTPUT+'.SNPS.FOUNDERS'


            TMP_merged_list = os.path.join(INTERMEDIATE_PATH, "merge_list")


            command = ' '.join(["echo", __AA__ + '.bed', __AA__ + '.bim', __AA__ + '.fam', ">", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(["echo", __HLA__ + '.bed', __HLA__ + '.bim', __HLA__ + '.fam', ">>", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join(["echo", __SNPS__ + '.bed', __SNPS__ + '.bim', __SNPS__ + '.fam', ">>", TMP_merged_list])
            # print(command)
            os.system(command)

            command = ' '.join([plink, "--merge-list", TMP_merged_list, "--make-bed", "--out", OUTPUT + '.MERGED'])
            # print(command)
            os.system(command)



            if not f_save_intermediates:

                os.system("rm " + __AA__+".bed")
                os.system("rm " + __AA__+".bim")
                os.system("rm " + __AA__+".fam")
                os.system("rm " + __AA__+".log")

                os.system("rm " + __SNPS__+".bed")
                os.system("rm " + __SNPS__+".bim")
                os.system("rm " + __SNPS__+".fam")
                os.system("rm " + __SNPS__+".log")

                os.system("rm " + __HLA__+".bed")
                os.system("rm " + __HLA__+".bim")
                os.system("rm " + __HLA__+".fam")
                os.system("rm " + __HLA__+".log")

                os.system("rm " + TMP_merged_list)

                if os.path.exists(__HLA__+".nosex"):
                    os.system("rm " + __HLA__ + ".nosex")
                if os.path.exists(__SNPS__+".nosex"):
                    os.system("rm " + __SNPS__ + ".nosex")
                if os.path.exists(__AA__+".nosex"):
                    os.system("rm " + __AA__ + ".nosex")

            index += 1


        if QC:

            """
            plink --bfile $OUTPUT.MERGED.FOUNDERS --freq --out $OUTPUT.MERGED.FOUNDERS.FRQ
            awk '{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > all.remove.snps
            awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order
        
            # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
            plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --exclude all.remove.snps --geno 0.5 --make-bed --out $OUTPUT
        
            # Calculate allele frequencies
            plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
            rm $SNP_DATA.FOUNDERS.*
            rm $OUTPUT.MERGED.FOUNDERS.*
            rm $OUTPUT.*.FOUNDERS.???
            rm allele.order
            rm all.remove.snps
            
            
            # 2020. 04. 29.
            This part is quite different from the one of HATK.
            
            """

            print("[{}] Performing quality control.".format(index))


            # Output of previous code block.
            __MERGED__ = OUTPUT + '.MERGED'
            TMP_allele_order = __MERGED__ + ".refallele"

            command = ' '.join(
                ["awk", '\'{if (NR > 1){if (($5 == "a" && $6 == "p") || ($6 == "a" && $5 == "p")){print $2 "\tp"}}}\'',
                 __MERGED__+'.bim', ">", TMP_allele_order])
            # print(command)
            os.system(command)

            command = ' '.join(
                [plink, "--make-bed", "--bfile", __MERGED__, "--a1-allele", TMP_allele_order,
                 "--out", OUTPUT]) # (2019. 01. 10.) Final output as just output prefix(`OUTPUT`)
            # print(command)
            os.system(command)


            if not f_save_intermediates:
                os.system("rm " + __MERGED__+".bed")
                os.system("rm " + __MERGED__+".bim")
                os.system("rm " + __MERGED__+".fam")
                os.system("rm " + __MERGED__+".log")
                os.system("rm " + __MERGED__+".refallele")


            index += 1


        # if PREPARE:
        #
        #     print("[{}] Preparing files for Beagle.".format(index))
        #
        #     """
        #     [Source from Buhm Han.]
        #
        #     awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
        #     plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
        #     awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
        #     cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
        #
        #     echo "[$i] Converting to beagle format.";  @ i++
        #     linkage2beagle pedigree=$OUTPUT.nopheno.ped data=$OUTPUT.dat beagle=$OUTPUT.bgl standard=true > $OUTPUT.bgl.log
        #
        #
        #     [Source from Yang.]
        #
        #     awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
        #     plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
        #     plink --bfile $OUTPUT --recode --transpose --out $OUTPUT
        #     # awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
        #     # cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped
        #
        #     echo "[$i] Converting to beagle format.";  @ i++
        #     beagle2vcf -fnmarker $OUTPUT.markers -fnfam $OUTPUT.fam -fngtype $OUTPUT.tped -fnout $OUTPUT.vcf
        #
        #     I will make this code block based on source given by Yang. for now.
        #
        #     """
        #
        #     # *.markers
        #     command = ' '.join(
        #         ["awk", '\'{print $2 " " $4 " " $5 " " $6}\'', OUTPUT + '.bim', ">", OUTPUT + '.markers'])
        #     print(command)
        #     os.system(command)
        #
        #     command = ' '.join(
        #         [plink, "--bfile", OUTPUT, "--keep-allele-order", "--recode", "--alleleACGT", "--out", OUTPUT])
        #     print(command)
        #     os.system(command)
        #
        #     # *.dat
        #     command = ' '.join(["awk", '\'{print "M " $2}\'', OUTPUT + '.map', ">", OUTPUT + '.dat'])
        #     print(command)
        #     os.system(command)
        #
        #     # *.nopheno.ped
        #     command = ' '.join(["cut -d ' ' -f1-5,7-", OUTPUT + '.ped', ">", OUTPUT + '.nopheno.ped'])
        #     print(command)
        #     os.system(command)
        #
        #     print("[{}] Converting to beagle format.".format(index))
        #
        #     command = ' '.join(
        #         [linkage2beagle, "pedigree=" + OUTPUT + '.nopheno.ped', "data=" + OUTPUT + '.dat',
        #          "beagle=" + OUTPUT + '.bgl', "standard=true", ">", OUTPUT + '.bgl.log'])
        #     print(command)
        #     os.system(command)
        #
        #
        #     index += 1
        #
        #
        #
        # if PHASE:
        #
        #     print("[{}] Phasing reference using Beagle (see progress in $OUTPUT.bgl.log).".format(index))
        #
        #     # Put this part postponed. (2017.11.29. by B. Han.)
        #     # Anyway introduced phasing by beagle. (2018. 7. 16.)
        #
        #     '''
        #     beagle unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log
        #
        #     '''
        #
        #     command = ' '.join(
        #         [beagle, "unphased=" + OUTPUT + '.bgl',
        #          "nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000",
        #          "log=" + OUTPUT + '.phasing', ">>", OUTPUT + '.bgl.log'])
        #     print(command)
        #     os.system(command)
        #
        #
        #     index += 1
        #
        #
        # if CLEANUP:
        #
        #     print("[{}] Removing unnecessary files.".format(index))
        #     '''
        #     rm $OUTPUT.nopheno.ped
        #     rm $OUTPUT.bgl.gprobs
        #     rm $OUTPUT.bgl.r2
        #     rm $OUTPUT.bgl
        #     rm $OUTPUT.ped
        #     rm $OUTPUT.map
        #     rm $OUTPUT.dat
        #     rm $OUTPUT.phasing.log
        #     '''
        #
        #     rm_tlist = ('.nopheno.ped', '.bgl.gprobs', '.bgl.r2', '.bgl', '.ped', '.map', '.dat')
        #
        #     for i in rm_tlist:
        #         print("rm " + OUTPUT + ".MERGED" + i) # (2019. 01. 10.) ".MERGED" should be removed.
        #         os.system("rm " + OUTPUT + ".MERGED" + i)
        #
        #     index += 1

        print("[{}] Making reference panel for HLA-AA,SNPS,HLA is Done!".format(index))



    return _OUT




# if __name__ == "__main__":
#
#     parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
#                                      description=textwrap.dedent('''\
#     #################################################################################################
#
#         MakeReference_v2.py
#
#         Generating markers based on HLA sequence information dictionary(generated by "IMGT2Seq").
#
#         (ex.)
#         : python3 MakeReference_v2.py
#             --chped
#             --hg
#             -o
#             --dict-AA
#             --dict-SNPS
#
#         HLA PED file should contain HLA alleles in the following (alphabetical) order:
#         HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1
#
#     #################################################################################################
#                                      '''),
#                                      # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
#                                      add_help=False)
#
#     parser._optionals.title = "OPTIONS"
#
#     parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
#
#     parser.add_argument("--variants", help="\nInput variants data file(.bed/.bim/.fam)\n\n")
#     parser.add_argument("--chped", help="\nHLA Type Data(.chped)\n\n", required=True)
#     parser.add_argument("--hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")
#     parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)
#
#     parser.add_argument("--dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", required=True)
#     parser.add_argument("--dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", required=True)
#
#     parser.add_argument("--save-intermediates", help="\nDon't remove intermediate files.\n\n", action='store_true')
#     parser.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')
#
#     # Beagle4.1
#     parser.add_argument("--phasing", help="\nPerform phasing with Beagle4.1.\n\n", action='store_true')
#     parser.add_argument("--mem", help="\nJava Memory requried for Bealge4.1. (ex. 2g)\n\n", default="2g")
#     parser.add_argument("--nthreads", help="\nThe number of threads to use in Bealge4.1. (ex. 2)\n\n", default=1, type=int)
#
#
#
#     ##### <for Test> #####
#
#     # 2019. 01. 10
#     # args = parser.parse_args(["-chped", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU_HLA.imgt370.4field.chped",
#     #                           "-variants", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HAPMAP_CEU",
#     #                           "-o", "/Users/wansun/Git_Projects/HATK/tests/_2_b_MarkerGenerator/20190110_bMarkerTest/HAPMAP_CEU_HLA.imgt370.hg18",
#     #                           "-hg", "18",
#     #                           "-dict-AA", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_AA.hg18.imgt370",
#     #                           "-dict-SNPS", "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/HATK/data/b_MarkerGenerator/HLA_DICTIONARY_SNPS.hg18.imgt370"
#     #                           ])
#
#     ##### <for Publication> #####
#
#     args = parser.parse_args()
#     print(args)
#
#
#     # Implementing Main Function.
#     MakeReference_v2(_CHPED=args.chped, _OUT=args.out, _hg=args.hg, _dictionary_AA=args.dict_AA,
#                      _dictionary_SNPS=args.dict_SNPS, _variants=args.variants, _mem=args.mem,
#                      f_save_intermediates=args.save_intermediates, f_phasing=args.phasing, _nthreads=args.nthreads)