# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
from pathlib import Path
from platform import platform


def MakeReference(_INPUT_DATA, _HLA_ped, _OUTPUT_Prefix,
                  _p_plink='./dependency/plink_mac' if not bool(re.search(pattern="Linux", string=platform())) else './dependency/plink_linux',
                  _p_beagle='./dependency/beagle.jar',
                  _p_linkage2beagle='./dependency/linkage2beagle.jar',
                  _dictionary_AA_map='./data/MakeReference_old/HLA_DICTIONARY_AA.map',
                  _dictionary_AA='./data/MakeReference_old/HLA_DICTIONARY_AA.txt',
                  _dictionary_SNPS_map='./data/MakeReference_old/HLA_DICTIONARY_SNPS.map',
                  _dictionary_SNPS='./data/MakeReference_old/HLA_DICTIONARY_SNPS.txt',
                  _previous_version=False,
                  _hg = "19",
                  _mem = "2000m"):



    ########## < Core Variables > ##########

    ### Module name

    std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
    print(std_MAIN_PROCESS_NAME + "Init.\n")


    ### Major paths for project

    # data
    p_data_MakeReference = './data/MakeReference_old' if _previous_version else './data/MakeReference'

    # src
    p_src_MakeReferece = "./src/MakeReference_old" if _previous_version else "./src/MakeReference"


    ########## <CHECK FOR DEPENDENCIES> ##########

    ### Other Software.

    if not os.path.exists(_p_plink):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(os.path.dirname(_p_plink)))
        sys.exit()
    if not os.path.exists(_p_beagle):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'Beagle 3' (http://faculty.washington.edu/browning/beagle/beagle.html#download) in '{0}'\n".format(os.path.dirname(_p_beagle)))
        sys.exit()
    if not os.path.exists(_p_linkage2beagle):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'linkage2beagle.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) (beagle.3.0.4/utility/linkage2beagle.jar) in '{0}'\n".format(os.path.dirname(_p_linkage2beagle)))
        sys.exit()

    ### Dictionary Information for HLA sequence

    if not os.path.exists(_dictionary_AA_map):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_AA.map' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_AA_map)))
        sys.exit()
    if not os.path.exists(_dictionary_AA):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_AA.txt' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_AA)))
        sys.exit()
    if not os.path.exists(_dictionary_SNPS_map):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_SNPS.map' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_SNPS_map)))
        sys.exit()
    if not os.path.exists(_dictionary_SNPS):
        print(std_MAIN_PROCESS_NAME + "Please Prepare 'HLA_DICTIONARY_SNPS.txt' (included with this package) in '{0}'\n".format(os.path.dirname(_dictionary_SNPS)))
        sys.exit()

    ### Source Code Scripts

    if not _previous_version:

        # New version with Python.

        if not os.path.exists(os.path.join(p_src_MakeReferece, "HLAtoSequences.py")):
            print(std_MAIN_PROCESS_NAME + "Error. 'HLAtoSequences.py' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()
        else:
            from src.MakeReference.HLAtoSequences import HLAtoSequences

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeVariants.py")):
            print(std_MAIN_PROCESS_NAME + "Error. 'encodeVariants.py' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()
        else:
            from src.MakeReference.encodeVariants import encodeVariants

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeHLA.py")):
            print(std_MAIN_PROCESS_NAME + "Error. 'encodeHLA.py' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()
        else:
            from src.MakeReference.encodeHLA import encodeHLA

    else:
        # Previous version with Perl.

        if not os.path.exists(os.path.join(p_src_MakeReferece, "HLAtoSequences.pl")):
            print(std_MAIN_PROCESS_NAME + "Error. 'HLAtoSequences.pl' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeVariants.pl")):
            print(std_MAIN_PROCESS_NAME + "Error. 'encodeVariants.pl' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeHLA.pl")):
            print(std_MAIN_PROCESS_NAME + "Error. 'encodeHLA.pl' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()



    ########## <Core Variables> ##########

    HLA_DATA = _HLA_ped


    # (2018. 5. 29) additional processing for intermediate path and output file prefix
    p = Path(os.path.abspath(_OUTPUT_Prefix))
    OUTPUT = _OUTPUT_Prefix
    INTERMEDIATE_PATH = os.path.join(*p.parts[:-1])


    # (2018. 5. 29) additional processing for intermediate path and output file prefix
    SNP_DATA = _INPUT_DATA
    _INPUT_DATA_prefix = Path(_INPUT_DATA).name

    SNP_DATA2 = os.path.join(INTERMEDIATE_PATH, _INPUT_DATA_prefix) # for each output directory

    print(std_MAIN_PROCESS_NAME + "Intermediate folder path is {0},\nOutput filename prefix is {1}\n\n".format(INTERMEDIATE_PATH, OUTPUT))

    # print(' '.join(["mkdir -p", INTERMEDIATE_PATH]))
    os.system(' '.join(["mkdir -p", INTERMEDIATE_PATH]))

    plink = ' '.join([_p_plink, "--noweb", "--silent"])
    beagle = ' '.join(["java", "-Xmx", _mem, "-jar", _p_beagle])
    linkage2beagle = ' '.join(["java", "-Xmx",_mem, "-jar", _p_linkage2beagle])



    ########## <Flags for Code Block> ##########

    ENCODE_AA = 1
    ENCODE_HLA = 1
    ENCODE_SNPS = 1
    EXTRACT_FOUNDERS = 1
    MERGE = 1
    QC = 1
    PREPARE = 1
    PHASE = 1
    CLEANUP = 0 # set to zero for time being



    ########## <Making Reference Panel> ##########

    print(std_MAIN_PROCESS_NAME + "Making Reference Panel for \"{0}\"".format(OUTPUT))

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

        print("\n[1] Generating amino acid sequences from HLA types.")

        if not _previous_version:
            HLAtoSequences(HLA_DATA, _dictionary_AA, "AA", _out=OUTPUT)
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "HLAtoSequences.pl"), HLA_DATA, _dictionary_AA, "AA", ">", OUTPUT+".AA.ped"])
            print(command)
            os.system(command)

        os.system(' '.join(["cp", _dictionary_AA_map, OUTPUT + '.AA.map']))

        print("\n[2] Encoding amino acids positions.")

        if not _previous_version:
            encodeVariants(OUTPUT + '.AA.ped', OUTPUT + '.AA.map', OUTPUT + '.AA.CODED') # previously "enCODED".
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "encodeVariants.pl"), OUTPUT+".AA.ped", OUTPUT+".AA.map", OUTPUT+".AA.CODED"])
            print(command)
            os.system(command)

        # command for checking output from encodeVariant.py(.pl)
        command = ' '.join([plink, "--file", OUTPUT+'.AA.CODED', "--missing-genotype 0", "--make-bed", "--out", OUTPUT+'.AA.TMP'])
        print(command)
        os.system(command)

        command = ' '.join(["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.AA.TMP.bim', "|", "grep -v {0}".format("INDEL" if not _previous_version else "INS"), "|", "cut -f2", ">", os.path.join(INTERMEDIATE_PATH, "to_remove")])

        """
        In previous framework originally created by Sherman Jia, Only insertions were dealt with as a marker "INS".
        In the new version of Framework, marker label is "INDEL".
        """

        print(command)
        os.system(command)

        command = ' '.join([plink, "--bfile", OUTPUT+'.AA.TMP', "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"), "--make-bed", "--out", OUTPUT+'.AA.CODED'])
        print(command)
        os.system(command)

        rm_tlist = (OUTPUT+'.AA.TMP*', os.path.join(INTERMEDIATE_PATH, "to_remove"), OUTPUT+'.AA.???')

        for i in rm_tlist:
            print(i)
            os.system("rm "+i)



    if ENCODE_HLA:

        print("\n[3] Encoding HLA alleles.")

        if not _previous_version:
            encodeHLA(HLA_DATA, OUTPUT, _hg)
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "encodeHLA.pl"), HLA_DATA, OUTPUT+".HLA.map", ">", OUTPUT+".HLA.ped"])
            print(command)
            os.system(command)


        command = ' '.join([plink, "--file", OUTPUT+'.HLA', "--make-bed", "--out", OUTPUT+'.HLA'])
        print(command)
        os.system(command)



    if ENCODE_SNPS:

        print("\n[4] Generating DNA sequences from HLA types.")

        if not _previous_version:
            HLAtoSequences(HLA_DATA, _dictionary_SNPS, "SNPS", OUTPUT)
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "HLAtoSequences.pl"), HLA_DATA, _dictionary_SNPS, "SNPS", ">", OUTPUT+".SNPS.ped"])
            print(command)
            os.system(command)

        command = ' '.join(["cp", _dictionary_SNPS_map, OUTPUT+'.SNPS.map'])
        print(command)
        os.system(command)

        print("\n[5] Encoding SNP positions.")

        if not _previous_version:
            encodeVariants(OUTPUT+'.SNPS.ped', OUTPUT+'.SNPS.map', OUTPUT+'.SNPS.CODED')
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "encodeVariants.pl"), OUTPUT+".SNPS.ped", OUTPUT+".SNPS.map", OUTPUT+".SNPS.CODED"])
            print(command)
            os.system(command)


        command = ' '.join([plink, "--file", OUTPUT+'.SNPS.CODED', "--missing-genotype 0", "--make-bed", "--out", OUTPUT+'.SNPS.TMP'])
        print(command)
        os.system(command)

        command = ' '.join(["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT +'.SNPS.TMP.bim', "|", "grep -v {0}".format("INDEL" if not _previous_version else "INS"), "|", "cut -f2", ">", os.path.join(INTERMEDIATE_PATH, "to_remove")])
        print(command)
        os.system(command)

        command = ' '.join([plink, "--bfile", OUTPUT+'.SNPS.TMP', "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"), "--make-bed", "--out", OUTPUT+'.SNPS.CODED'])
        print(command)
        os.system(command)


        rm_tlist = (OUTPUT+'.SNPS.TMP*', os.path.join(INTERMEDIATE_PATH, 'to_remove'), OUTPUT+'.SNPS.???')

        for i in rm_tlist:
            print(i)
            os.system("rm "+i)



    if EXTRACT_FOUNDERS:

        print("\n[5] Encoding SNP positions.")

        """
        if ($EXTRACT_FOUNDERS) then
            echo "[$i] Extracting founders."; @ i++
            # founder의 정의가 이게 맞는건 아니겠지만, 아래 plink명령어를 거치고 나오는 founder라는 애들은 모두 엄마,아빠 ID정보가 없는 애들임.
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

        command = ' '.join([plink, "--bfile", SNP_DATA, "--filter-founders", "--mind 0.3", "--alleleACGT", "--make-bed", "--out", os.path.join(INTERMEDIATE_PATH, _INPUT_DATA_prefix+'.FOUNDERS')])
        print(command)
        os.system(command)

        # SNP_DATA2 = os.path.join(INTERMEDIATE_PATH, _INPUT_DATA_prefix)

        # Initial QC on Reference SNP panel
        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--hardy", "--out", SNP_DATA2+'.FOUNDERS.hardy'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--freq", "--out", SNP_DATA2+'.FOUNDERS.freq'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--missing", "--out", SNP_DATA2+'.FOUNDERS.missing'])
        print(command)
        os.system(command)

        command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.hardy.hwe', "|", "awk", "' $9 < 0.000001 { print $2 }'", "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.hardy")])
        print(command)
        os.system(command)
        command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.freq.frq', "|", "awk", "' $5 < 0.01 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.freq")])
        print(command)
        os.system(command)
        command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.missing.lmiss', "|", "awk", "' $5 > 0.05 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.missing")])
        print(command)
        os.system(command)
        command = ' '.join(["cat", os.path.join(INTERMEDIATE_PATH, "remove.snps.*"), "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "all.remove.snps")])
        print(command)
        os.system(command)


        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--allow-no-sex", "--exclude", os.path.join(INTERMEDIATE_PATH, "all.remove.snps"), "--make-bed", "--out", SNP_DATA2+'.FOUNDERS.QC'])
        print(command)
        os.system(command)

        # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

        command = ' '.join([plink, "--bfile", OUTPUT+'.HLA', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.HLA.FOUNDERS'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", OUTPUT+'.SNPS.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.SNPS.FOUNDERS'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", OUTPUT+'.AA.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.AA.FOUNDERS'])
        print(command)
        os.system(command)

        command = ' '.join(["rm", os.path.join(INTERMEDIATE_PATH, "remove.snps.*")])
        print(command)
        os.system(command)



    if MERGE:

        print("\n[6] Merging SNP, HLA, and amino acid datasets.")

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

        TMP_merged_list = os.path.join(INTERMEDIATE_PATH, "merge_list")

        command = ' '.join(["echo", OUTPUT+'.HLA.FOUNDERS.bed', OUTPUT+'.HLA.FOUNDERS.bim', OUTPUT+'.HLA.FOUNDERS.fam', ">", TMP_merged_list])
        print(command)
        os.system(command)

        command = ' '.join(["echo", OUTPUT+'.AA.FOUNDERS.bed', OUTPUT+'.AA.FOUNDERS.bim', OUTPUT+'.AA.FOUNDERS.fam', ">>", TMP_merged_list])
        print(command)
        os.system(command)

        command = ' '.join(["echo", OUTPUT+'.SNPS.FOUNDERS.bed', OUTPUT+'.SNPS.FOUNDERS.bim', OUTPUT+'.SNPS.FOUNDERS.fam', ">>", TMP_merged_list])
        print(command)
        os.system(command)

        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS.QC', "--merge-list", TMP_merged_list, "--make-bed", "--out", OUTPUT+'.MERGED.FOUNDERS'])
        print(command)
        os.system(command)

        rm_tlist = (OUTPUT+'.HLA.???', OUTPUT+'.AA.CODED.???', OUTPUT+'.SNPS.CODED.???', TMP_merged_list)

        for i in rm_tlist:
            print(i)
            os.system("rm "+i)



    if QC:

        print("\n[7] Performing quality control.")

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

        TMP_allele_order = os.path.join(INTERMEDIATE_PATH, "allele.order")
        TMP_all_remove_snps = os.path.join(INTERMEDIATE_PATH, "all.remove.snps")

        command = ' '.join([plink, "--bfile", OUTPUT+'.MERGED.FOUNDERS', "--freq", "--out", OUTPUT+'.MERGED.FOUNDERS.FRQ'])
        print(command)
        os.system(command)
        command = ' '.join(["awk", "'{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}'", OUTPUT+'.MERGED.FOUNDERS.FRQ.frq', ">", TMP_all_remove_snps])
        print(command)
        os.system(command)
        command = ' '.join(["awk", '\'{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}\'', OUTPUT+'.MERGED.FOUNDERS.FRQ.frq', ">", TMP_allele_order])
        print(command)
        os.system(command)

        # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
        command = ' '.join([plink, "--bfile", OUTPUT+'.MERGED.FOUNDERS', "--reference-allele", TMP_allele_order, "--exclude", TMP_all_remove_snps, "--geno 0.5", "--make-bed", "--out", OUTPUT])
        print(command)
        os.system(command)

        # Calculate allele frequencies
        command = ' '.join([plink, "--bfile", OUTPUT, "--keep-allele-order", "--freq", "--out", OUTPUT+'.FRQ'])
        print(command)
        os.system(command)

        rm_tlist = (SNP_DATA2+'.FOUNDERS.*', OUTPUT+'.MERGED.FOUNDERS.*', OUTPUT+'.*.FOUNDERS.???', TMP_allele_order, TMP_all_remove_snps)

        for i in rm_tlist:
            print(i)
            os.system("rm "+i)



    if PREPARE:

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

        print("\n[8] Preparing files for Beagle.")

        command = ' '.join(["awk", '\'{print $2 " " $4 " " $5 " " $6}\'', OUTPUT+'.bim', ">", OUTPUT+'.markers'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", OUTPUT, "--keep-allele-order", "--recode", "--alleleACGT", "--out", OUTPUT])
        print(command)
        os.system(command)
        command = ' '.join(["awk", '\'{print "M " $2}\'', OUTPUT+'.map', ">", OUTPUT+'.dat'])
        print(command)
        os.system(command)
        command = ' '.join(["cut -d ' ' -f1-5,7-", OUTPUT+'.ped', ">", OUTPUT+'.nopheno.ped'])
        print(command)
        os.system(command)

        print("\n[9] Converting to beagle format.")

        command = ' '.join([linkage2beagle, "pedigree="+OUTPUT+'.nopheno.ped', "data="+OUTPUT+'.dat', "beagle="+OUTPUT+'.bgl', "standard=true", ">", OUTPUT+'.bgl.log'])
        print(command)
        os.system(command)



    if PHASE:

        # Put this part postponed. (2017.11.29. by B. Han.)
        # Anyway introduced phasing by beagle. (2018. 7. 16.)

        '''
        beagle unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log

        '''
        print("\n[10] Phasing reference using Beagle (see progress in $OUTPUT.bgl.log).")

        command= ' '.join([beagle, "unphased="+OUTPUT+'.bgl', "nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000", "log="+OUTPUT+'.phasing', ">>", OUTPUT+'.bgl.log'])
        print(command)
        os.system(command)



    if CLEANUP:

        print("\n[11] Removing unnecessary files.")
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


    print("\n[12] Done!")




    return 0


if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        MakeReference.py

        This script helps prepare a reference dataset for HLA imputation

        Usage(1)
        : python3 MakeReference.py --previous-version -i ./data/MakeReference_old/HAPMAP_CEU
            -ped ./data/MakeReference_old/HAPMAP_CEU_HLA.ped -hg 18 -o ./Trial_HAPMAP_CEU

        Usage(2)
        : python3 MakeReference.py -i ./data/MakeReference/HAPMAP_CEU
            -ped ./data/MakeReference/HAPMAP_CEU_HLA.4field.ped -hg 18 -o ./Trial_HAPMAP_CEU
            -dict-AA ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt
            -dict-AA-map ./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map
            -dict-SNPS ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt
            -dict-SNPS-map ./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map

        HLA PED file should contain HLA alleles in the following (alphabetical) order:
        HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1

    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-i", help="\nInput Data file(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("-ped", help="\nHLA Type Data(.ped)\n\n", required=True)
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")
    parser.add_argument("-mem", help="\nMemory requried for beagle\n\n", default="12g")
    parser.add_argument("-o", help="\nOutput file prefix\n\n")

    parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as old version.\n\n",
                        action='store_true')

    hla_dict = parser.add_argument_group(title='HLA_DICTIONARY',
                                         description='- Arguments to specify HLA_DICTIONARY Information to New version of MakeReference(2.0)\n'
                                                     '- If you\'re going to use previous version of MakeReference, then Don\'t care about these options.')

    hla_dict.add_argument("-dict-AA", help="\nInput HLA Dictionary file for AA Information.\n\n", default="Not_given")
    hla_dict.add_argument("-dict-AA-map", help="\nInput HLA Dictionary .map file for AA Information.\n\n", default="Not_given")
    hla_dict.add_argument("-dict-SNPS", help="\nInput HLA Dictionary file for SNPS Information.\n\n", default="Not_given")
    hla_dict.add_argument("-dict-SNPS-map", help="\nInput HLA Dictionary .map file for SNPS Information\n\n", default="Not_given")




    ##### <for Publication> #####

    args = parser.parse_args()


    ##### <for Test> #####

    # ==========< New version >==========

    # # v370, hg18 (perfectly what Sherman dealt with.)
    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU",
    #                           "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-hg", "18",
    #                           "-o", "./MAKEREFERENCE/MAKEREFERENCE_PYTHON",
    #                           "-dict-AA", "./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-dict-AA-map", "./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map",
    #                           "-dict-SNPS", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-dict-SNPS-map", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map",
    #                           ])

    # # v3320
    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU",
    #                           "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.ped",
    #                           "-o", "OLD_VERSION_TEST/PREV_VERSION_TEST",
    #                           "-dict-AA", "./data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.txt",
    #                           "-dict-AA-map", "./data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.map",
    #                           "-dict-SNPS", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg19.imgt3320.txt",
    #                           "-dict-SNPS-map", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg19.imgt3320.map",
    #"-m","2000m",
    #                           ])

    # ==========< Perfectly Old version >==========

    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU", "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.ped", "-o", "PREV_VERSION_TEST", "--previous-version"]) # 완전히 구버젼으로 작동하게 하고 싶을때.

    # # Intermediate Path
    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU", "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.ped", "-o", "OLD_VERSION_TEST/PREV_VERSION_TEST", "--previous-version"]) # 완전히 구버젼으로 작동하게 하고 싶을때.


    ##### Additional Argument processing

    t_dict_AA = ""
    t_dict_AA_map = ""
    t_dict_SNPS = ""
    t_dict_SNPS_map = ""
    t_mem = ""


    if (args.dict_AA != "Not_given" and args.dict_AA_map != "Not_given" and args.dict_SNPS != "Not_given" and args.dict_SNPS_map != "Not_given"):

        # When all HLA DICTIONARY information is given properly,

        t_dict_AA = args.dict_AA
        t_dict_AA_map = args.dict_AA_map
        t_dict_SNPS = args.dict_SNPS
        t_dict_SNPS_map = args.dict_SNPS_map

    elif (args.dict_AA == "Not_given" and args.dict_AA_map == "Not_given" and args.dict_SNPS == "Not_given" and args.dict_SNPS_map == "Not_given"):

        # No values are given to HLA DICTIONARY related options.

        if not args.previous_version:
            # Abort
            print("\n[Error]: None of HLA DICTIONARY files are given. Please check them all again.")
            print('{"-dict-AA", "-dict-AA-map", "-dict-SNPS", "-dict-SNPS-map"}\n')
            sys.exit()

        else:
            # In case of old version of MakeReference, then use default dataset.
            t_dict_AA = './data/MakeReference_old/HLA_DICTIONARY_AA.txt'
            t_dict_AA_map = './data/MakeReference_old/HLA_DICTIONARY_AA.map'
            t_dict_SNPS = './data/MakeReference_old/HLA_DICTIONARY_SNPS.txt'
            t_dict_SNPS_map = './data/MakeReference_old/HLA_DICTIONARY_SNPS.map'

    else:
        # Abort
        print("\n[Error]: Not all of HLA DICTIONARY files are given. Please check them all again.")
        print('{"-dict-AA", "-dict-AA-map", "-dict-SNPS", "-dict-SNPS-map"}\n')
        sys.exit()


    # (2018. 7. 16.)
    if args.previous_version:

        """
        If the option `--previous-version` is given, then the value of `-hg` option will be fixed to "18".
        I took this measure because previous version of the framework only covers hg 18.
        """

        args.hg = "18"
        args.mem= "2000m"

        print("\n[Warning]: Only hg18 is available for previous version of MakeReference. Human Genome version will be set to hg 18 by force.\n")




    print(args)


    # Implementing Main Function.
    MakeReference(_INPUT_DATA=args.i, _HLA_ped=args.ped, _OUTPUT_Prefix=args.o, _previous_version=args.previous_version, _hg=args.hg,
                  _dictionary_AA=t_dict_AA, _dictionary_AA_map=t_dict_AA_map,
                  _dictionary_SNPS=t_dict_SNPS, _dictionary_SNPS_map=t_dict_SNPS_map,
                  _mem=args.mem)
