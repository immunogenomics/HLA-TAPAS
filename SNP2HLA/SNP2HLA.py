# -*- coding: utf-8 -*-

import os, sys, re
import subprocess
from shutil import which


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

def exec_cmd(_command):
    subprocess.run(_command, shell=True, check=True)



def SNP2HLA(_target, _reference, _out, _mem="2g", _tolerated_diff=.15, _p_dependency="./dependency",
            _b_nthreads=1, _b_niterations=5):


    ### [0] Input data

    # Target data check.
    if not os.path.exists(_target+'.bed'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Target SNP data can't be found('{}'). Please check '--target/-t' argument again.".format(_target + '.bed'))
        sys.exit()
    if not os.path.exists(_target+'.bim'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Target SNP data can't be found('{}'). Please check '--target/-t' argument again.".format(_target + '.bim'))
        sys.exit()
    if not os.path.exists(_target+'.fam'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Target SNP data can't be found('{}'). Please check '--target/-t' argument again.".format(_target + '.fam'))
        sys.exit()

    # Reference data check.
    # if not os.path.exists(_reference+'.bed'):
    #     print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.bed'))
    #     sys.exit()
    if not os.path.exists(_reference+'.bim'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.bim'))
        sys.exit()
    # if not os.path.exists(_reference+'.fam'):
    #     print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.fam'))
    #     sys.exit()
    if not os.path.exists(_reference+'.FRQ.frq'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.FRQ.frq'))
        sys.exit()
    if not os.path.exists(_reference+'.markers'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.bglv4.markers'))
        sys.exit()
    if not os.path.exists(_reference+'.bgl.phased.vcf.gz'):
        print(std_ERROR_MAIN_PROCESS_NAME + "One of Reference data can't be found('{}'). Please check '--reference/-ref' argument again.".format(_reference + '.bglv4.bgl.phased.vcf.gz'))
        sys.exit()




    ### [1] Major Path Variable and Dependent software

    if not os.path.exists(_p_dependency):
        print(std_ERROR_MAIN_PROCESS_NAME + "Folder for dependent software('{}') can't be found. Please check '--dependency' argument again.".format(_p_dependency))
        sys.exit()


    # dependent software.
    _plink = os.path.join(_p_dependency, "plink") #plink v1.9
    _beagle = os.path.join(_p_dependency, "beagle.jar")   # Beagle(v4.1)
    _linkage2beagle = os.path.join(_p_dependency, "linkage2beagle.jar")
    _beagle2vcf = os.path.join(_p_dependency, "beagle2vcf.jar")
    #_beagle2linkage = os.path.join(_dependency, "beagle2linkage.jar")
    _vcf2gprobs = os.path.join(_p_dependency, "vcf2gprobs.jar")
    _perl = which('perl')
    _merge_table = os.path.join("SNP2HLA/src/merge_tables.pl")
    # _csh = which('csh') if bool(which('csh')) else which('tcsh')
    # _parse_dosage = os.path.join("SNP2HLA/src/ParseDosage.csh")


    if not os.path.exists(_plink):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{0}'\n".format(_p_dependency))
        sys.exit()
    if not os.path.exists(_beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'Beagle 4.1' (https://faculty.washington.edu/browning/beagle/b4_1.html#download) in '{0}'\n".format(_p_dependency))
        sys.exit()
    if not os.path.exists(_linkage2beagle):
       print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'linkage2beagle.jar' in '{}'.".format(_p_dependency))
       sys.exit()
    if not os.path.exists(_beagle2vcf):
       print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'beagle2vcf.jar' in '{}'.".format(_p_dependency))
       sys.exit()
    #if not os.path.exists(_beagle2linkage):
    #    print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'beagle2linkage.jar' in 'dependency/' folder.")
    #    sys.exit()
    if not bool(_perl):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please install/prepare 'Perl' in your system.")
        sys.exit()
    # if not bool(_csh):
    #     print(std_ERROR_MAIN_PROCESS_NAME + "Please install/prepare 'csh/tcsh' in your system.")
    #     sys.exit()
    if not os.path.exists(_merge_table):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'merge_tables.pl' in 'SNP2HLA/src/' folder.")
        sys.exit()
    # if not os.path.exists(_parse_dosage):
    #     print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'ParseDosage.csh' in 'SNP2HLA/src/' folder.")
    #     sys.exit()



    ### [2] Memory representation check.
    p_Mb = re.compile(r'\d+m')
    p_Gb = re.compile(r'\d+[gG]')

    if not (bool(p_Mb.match(_mem)) or bool(p_Gb.match(_mem))):
        print(std_ERROR_MAIN_PROCESS_NAME + "Given Java memory value('{}') has bizzare representation. Please check '--mem' argument again.".format(_mem))
        sys.exit()



    ### [3] Intermediate path.

    OUTPUT = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        # If `os.path.dirname(OUTPUT)` doesn't exist, then it means the output of MakeReference should be genrated in current directory.
        INTERMEDIATE_PATH = "./"


    JAVATMP = _out+".javatmpdir"
    exec_cmd("mkdir -p " + JAVATMP)



    ### [4] Setting commands

    PLINK = ' '.join([_plink, "--silent", "--allow-no-sex"]) # "--noweb" won't be included because it is Plink1.9
    BEAGLE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _beagle])
    LINKAGE2BEAGLE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _linkage2beagle])
    BEAGLE2VCF = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _beagle2vcf])
    #BEAGLE2LINKAGE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _beagle2linkage])

    MERGE = ' '.join(["perl", _merge_table])
    #MERGE = _merge_table
    # PARSEDOSAGE = _parse_dosage



    ########## <Flags for Code Block> ##########
    EXTRACT_MHC = 1
    FLIP = 1
    IMPUTE = 1


    print("SNP2HLA: Performing HLA imputation for dataset {}".format(_target))
    print("- Java memory = {}b".format(_mem))


    index= 1
    __MHC__ = _out+".MHC"


    if EXTRACT_MHC:

        print("[{}] Extracting SNPs from the MHC.".format(index)); index += 1
        #MAF >1% as imputation threshold
        command = ' '.join([PLINK, "--bfile", _target, "--chr 6", "--from-mb 28 --to-mb 34", "--maf 0.01", "--make-bed", "--out", __MHC__])
        # print(command)
        exec_cmd(command)


    if FLIP:

        print("[{}] Performing SNP quality control.".format(index)); index += 1

        ### Identifying non-A/T non-C/G SNPs to flip
        command = ' '.join(["echo", "SNP 	POS	A1	A2", ">", OUTPUT+".tmp1"])
        # print(command)
        exec_cmd(command)
        command = ' '.join(["cut", "-f2,4-", __MHC__+".bim", ">>", OUTPUT+".tmp1"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["echo", "SNP    POSR    A1R A2R", ">", OUTPUT+".tmp2"])
        # print(command)
        exec_cmd(command)
        command = ' '.join(["cut", "-f2,4-", _reference + ".bim", ">>", OUTPUT + ".tmp2"])
        # print(command)
        exec_cmd(command)

        command = ' '.join([MERGE, OUTPUT+".tmp2", OUTPUT+".tmp1", "SNP", "|", "grep -v -w NA", ">", OUTPUT+".SNPS.alleles"])
        # print(command)
        exec_cmd(command)



        ### < Major flip 1 > ###

        command = ' '.join(["awk", "'{if ($3 != $6 && $3 != $7){print $1}}'", OUTPUT+".SNPS.alleles", ">", OUTPUT+".SNPS.toflip1"])
        # print(command)
        exec_cmd(command)

        command = ' '.join([PLINK, "--bfile", __MHC__, "--flip", OUTPUT+".SNPS.toflip1", "--make-bed", "--out", __MHC__+".FLP"])
        # print(command)
        exec_cmd(command)

        ## Calculating allele freqeuncy
        command = ' '.join([PLINK, "--bfile", __MHC__+".FLP", "--freq", "--out", __MHC__+".FLP.FRQ"])
        # print(command)
        exec_cmd(command)


        command = ' '.join(["sed 's/A1/A1I/g'", __MHC__+".FLP.FRQ.frq", "|", "sed 's/A2/A2I/g'", "|", "sed 's/MAF/MAF_I/g'", ">", OUTPUT+".tmp"])
        # print(command)
        exec_cmd(command)



        command = ' '.join(["mv", OUTPUT+".tmp", __MHC__+".FLP.FRQ"])
        # print(command)
        exec_cmd(command)

        command = ' '.join([MERGE, _reference + ".FRQ.frq", __MHC__ + ".FLP.FRQ.frq", "SNP", "|", "grep -v -w NA", ">", OUTPUT + ".SNPS.frq"])
        # print(command)
        exec_cmd(command)




        ### < Major flip 2 > ### (*.parsed file)
        command = ' '.join(["sed 's/ /\t/g'", OUTPUT+".SNPS.frq", "|",
                            'awk \'{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}\'',
                            ">", OUTPUT+".SNPS.frq.parsed"])
        # print(command)
        exec_cmd(command)



        ### < Major flip 3 > ###
        # Finding A/T and C/G SNPs
        command = ' '.join(['awk \'{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}\'',
                            OUTPUT+".SNPS.frq.parsed", ">", OUTPUT+".SNPS.ATCG.frq"])
        # print(command)
        exec_cmd(command)



        ### < Major flip 4 > ###

        # Identifying A/T and C/G SNPs to flip or remove
        command = ' '.join(["awk '{if ($10 < $9 && $10 < .15){print $1}}'", OUTPUT+".SNPS.ATCG.frq", ">", OUTPUT+".SNPS.toflip2"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["awk '{if ($4 > 0.4){print $1}}'", OUTPUT+".SNPS.ATCG.frq", ">", OUTPUT+".SNPS.toremove"])
        # print(command)
        exec_cmd(command)


        ## Identifying non A/T and non C/G SNPs to remove
        command = ' '.join(['awk \'{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > \'%f\'){print $1}}}\''%(_tolerated_diff),
                            OUTPUT+".SNPS.frq.parsed", ">>", OUTPUT+".SNPS.toremove"])
        # print(command)
        exec_cmd(command)


        command = ' '.join(['awk \'{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}\'',
                            OUTPUT+".SNPS.frq.parsed", ">>", OUTPUT+".SNPS.toremove"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(['awk \'{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}\'',
                            OUTPUT+".SNPS.frq.parsed", ">>", OUTPUT+".SNPS.toremove"])
        # print(command)
        exec_cmd(command)
        
        command = ' '.join(["sort", OUTPUT+".SNPS.toremove", "|","uniq > temp" ])
        # print(command)
        exec_cmd(command)
 
        command = ' '.join(['mv temp',OUTPUT+".SNPS.toremove" ] )
        # print(command)
        exec_cmd(command)

        ## Making QCd SNP file
        command = ' '.join([PLINK, "--bfile", __MHC__+".FLP", "--geno 0.2", "--exclude", OUTPUT+".SNPS.toremove", "--flip", OUTPUT+".SNPS.toflip2", "--make-bed", "--out", __MHC__+".QC"]) ###
        # print(command)
        exec_cmd(command)

        command = ' '.join([PLINK, "--bfile", __MHC__+".QC", "--freq", "--out", __MHC__+".QC.FRQ"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["sed 's/A1/A1I/g'", __MHC__+".QC.FRQ.frq", "|", "sed 's/A2/A2I/g'", "|", "sed 's/MAF/MAF_I/g'", ">", OUTPUT+".tmp"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["mv", OUTPUT+".tmp", __MHC__+".QC.FRQ.frq"])
        # print(command)
        exec_cmd(command)

        command = ' '.join([MERGE, _reference + ".FRQ.frq", __MHC__ + ".QC.FRQ.frq", "SNP", "|", "grep -v -w NA", ">", OUTPUT + ".SNPS.QC.frq"])
        # print(command)
        exec_cmd(command)


        command = ' '.join(["cut -f2", OUTPUT+".SNPS.QC.frq", "|", "awk '{if (NR > 1){print $1}}'", ">", OUTPUT+".SNPS.toinclude"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(['echo "SNP     POS    A1    A2"', ">", OUTPUT+".tmp1"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["cut -f2,4-", __MHC__+".QC.bim", ">>", OUTPUT+".tmp1"]) ######
        # print(command)
        exec_cmd(command)

        command = ' '.join([MERGE, OUTPUT+".tmp2", OUTPUT+".tmp1", "SNP", "|", 'awk \'{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}\'',
                            ">", __MHC__+".QC.bim"])
        # print(command)
        exec_cmd(command)



        ## Remove temporary files.
        exec_cmd(' '.join(["rm ", OUTPUT + ".tmp1"]))
        exec_cmd(' '.join(["rm ", OUTPUT + ".tmp2"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".FLP.*"]))
        # exec_cmd(' '.join(["rm ", __MHC__+".QC.ped"]))
        # exec_cmd(' '.join(["rm ", __MHC__+".QC.map"]))

        exec_cmd(' '.join(["rm ", __MHC__ + ".bed"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".bim"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".fam"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".log"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".filtered.bed"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".filtered.bim"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".filtered.fam"]))
        exec_cmd(' '.join(["rm ", __MHC__ + ".filtered.log"]))



        ### < Making *.vcf file for imputation (Beagle v4.x.x.) > ###

        """
        # Extracting SNPs and recoding QC'd file as vcf
        plink --bfile $MHC.QC --extract $OUTPUT.SNPS.toinclude --make-bed --out $MHC.QC.reorder
        plink --bfile $MHC.QC.reorder --recode vcf-iid --a1-allele $REFERENCE.markers 4 1 --out $MHC.QC
        
        """

        # Extracting SNPs and recoding QC'd file as vcf
        command = ' '.join([PLINK, "--bfile", __MHC__+".QC", "--extract", OUTPUT+".SNPS.toinclude", "--make-bed --out", __MHC__+".QC.reorder" ])
        # print(command)
        exec_cmd(command)


        exec_cmd(' '.join(["rm ", OUTPUT + ".SNPS.*"]))
        exec_cmd(' '.join(["rm ", __MHC__+".QC.bed"]))
        exec_cmd(' '.join(["rm ", __MHC__+".QC.bim"]))
        exec_cmd(' '.join(["rm ", __MHC__+".QC.fam"]))
        exec_cmd(' '.join(["rm ", __MHC__+".QC.log"]))
        exec_cmd(' '.join(["rm ", __MHC__+".QC.FRQ.frq"]))
        exec_cmd(' '.join(["rm ", __MHC__+".QC.FRQ.log"]))



        ## PLINK to Beagle format
        command = ' '.join(["awk", '\'{print $2" "$4" "$6" "$5}\'', __MHC__+".QC.reorder.bim", ">", __MHC__+".QC.markers"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["awk", '\'{print $2" "$5}\'', __MHC__+".QC.reorder.bim", ">", __MHC__+".QC.reorder.a1_allele"])
        # print(command)
        exec_cmd(command)

        """
        Note that '{print $2" "$4" "$6" "$5}' is used not '{print $2" "$4" "$5" "$6}'
        """

        command = ' '.join(
            [PLINK, "--bfile", __MHC__+".QC.reorder",
             "--keep-allele-order", "--recode", "--alleleACGT", "--out", __MHC__+".QC.reorder",
             "--a1-allele", __MHC__+".QC.reorder.a1_allele"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["awk", '\'{print "M " $2}\'', __MHC__+".QC.reorder.map", ">", __MHC__+".QC.reorder.dat"])
        # print(command)
        exec_cmd(command)

        command = ' '.join(["cut -d ' ' -f1-5,7-", __MHC__+".QC.reorder.ped", ">", __MHC__+".QC.reorder.nopheno.ped"])
        # print(command)
        exec_cmd(command)


        print("[{}] Converting PLINK to BEAGLE format.".format(index)); index += 1

        command = ' '.join(
            [LINKAGE2BEAGLE, "pedigree={}".format(__MHC__+".QC.reorder.nopheno.ped"), "data={}".format(__MHC__+".QC.reorder.dat"),
             "beagle={}".format(__MHC__+".QC.bgl"), "standard=true"])
        # print(command)

        try:
            # exec_cmd(command)
            f_log = open(__MHC__+".QC.bgl.log", 'w')
            subprocess.run(re.split(r'\s+', command), check=True, stdout=f_log, stderr=f_log)

        except subprocess.CalledProcessError:
            # Fail
            print(std_ERROR_MAIN_PROCESS_NAME + "'linkage2beagle.jar' failed. See log file('{}').".format(__MHC__+".QC.bgl.log"))
            sys.exit()
        else:
            # Succeed
            f_log.close()
            exec_cmd(' '.join(["rm ", __MHC__+".QC.reorder.*"]))
            exec_cmd(' '.join(["rm ", __MHC__+".QC.bgl.log"]))



        ## Beagle to VCF format
        print("[{}] Converting BEAGLE to VCF format.".format(index)); index += 1

        command = ' '.join([BEAGLE2VCF, '6', __MHC__+".QC.markers", __MHC__+".QC.bgl", '0', '>', __MHC__+".QC.bgl.vcf"]) # => Converted Target file to impute.
        # print(command)
        exec_cmd(command)

        exec_cmd(' '.join(["rm ", __MHC__+".QC.bgl"]))
        exec_cmd(' '.join(["rm ", __MHC__+".QC.markers"]))



    if IMPUTE:

        """
        if ($#argv >= 8) then
            beagle ref=$REFERENCE.bgl.vcf.gz gt=$MHC.QC.vcf impute=true gprobs=true nthreads=$THREAD chrom=6 niterations=$ITER lowmem=true out=$OUTPUT.bgl map=$MAP
	    else
            beagle ref=$REFERENCE.bgl.vcf.gz gt=$MHC.QC.vcf impute=true gprobs=true nthreads=$THREAD chrom=6 niterations=$ITER lowmem=true out=$OUTPUT.bgl
        """

        print("[{}] Performing HLA imputation.".format(index)); index += 1

        ## new beagle (>v4), assuming 4 threads and 10 interations
        command = ' '.join([BEAGLE,
                          "gt={}".format(__MHC__+".QC.bgl.vcf"),
                          'ref={}'.format(_reference + ".bgl.phased.vcf.gz"),
                          'impute=true',
                          'gprobs=true',
                          'nthreads={}'.format(_b_nthreads),
                          'chrom=6',
                          'niterations={}'.format(_b_niterations),
                          'lowmem=true',
                          'out={}'.format(OUTPUT +".bgl.phased")])

        # print(command)

        try:
            f_log = open(OUTPUT +".bgl.phased.vcf.log", 'w')
            subprocess.run(re.split(r'\s+', command), check=True, stdout=f_log, stderr=f_log)

        except subprocess.CalledProcessError:
            # Fail
            print(std_ERROR_MAIN_PROCESS_NAME + "Imputation failed. See log file('{}').".format(OUTPUT +".bgl.phased.vcf.log"))
            sys.exit()
        else:
            # Succeed
            f_log.close()
            exec_cmd(' '.join(["rm ", __MHC__+".QC.bgl.vcf"]))

            __IMPUTED__ = OUTPUT + ".bgl.phased.vcf.gz"
            return __IMPUTED__
        finally:
            # remove redundant log file.
            exec_cmd("rm {}".format(OUTPUT +".bgl.phased.log"))
            exec_cmd('rm -rf {}'.format(JAVATMP))

        # """
        # (1) Imputation result in *.vcf.gz file
        # (2) Imputation result in *.{bed,bim,fam} files (*.vcf.gz => *.{bed,bim,fam})
        # (2) Dosage file (*.gprobs => *.dosage)
        # """
        #
        #
        # # (2) Imputation result in *.{bed,bim,fam} files (*.vcf.gz => *.{bed,bim,fam})
        # command = ' '.join([PLINK, "--make-bed", "--vcf", __IMPUTED__, "--a1-allele {} 4 1".format(_reference + ".markers"), "--out", OUTPUT])
        # #print(command)
        # #exec_cmd(command)
        #
        #
        # # (3) Dosage file
        # command = ' '.join(["gunzip -c", __IMPUTED__, "|", "cat", "|", "java -jar {} > {}".format(_vcf2gprobs, OUTPUT+".bgl.gprobs")])
        # #print(command)
        # #exec_cmd(command)
        #
        # __gprobs__ = OUTPUT+".bgl.gprobs"
        #
        #
        # command = ' '.join(["tail -n +2 {}".format(__gprobs__), "|",
        #                     PARSEDOSAGE, "- > {}".format(OUTPUT+".dosage")])
        # #print(command)
        # #exec_cmd(command)
        #
        #
        #
        # exec_cmd(' '.join(["rm ", __MHC__+".QC.vcf.gz"]))
        # exec_cmd(' '.join(["rm -rf", JAVATMP]))

    print("Done\n")



# if __name__ == "__main__" :
#
#     parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
#                                      add_help=False,
#                                      prog='SNP2HLA',
#                                      description=textwrap.dedent('''\
#     #################################################################################################
#
#         < SNP2HLA.py >
#
#         SNP2HLA: Imputation of HLA amino acids and classical alleles from SNP genotypes
#
#         Author: Sherman Jia (xiaomingjia@gmail.com)
#                 + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
#                 + Extensive modifications by Phil Stuart (pstuart@umich.edu) on 4/19/16 to allow use of Beagle 4.1:
#                     verified to work with 22Feb16, 15Apr16, and 03May16 versions of Beagle.
#                 + Small modifications by Yang Luo (yangluo@broadinstitute.org): 09/30/16: verified working with Bealge 4.1 27Jun16.
#                 + Recoded to Python and updated by Wanson Choi(wansonchoi@snu.ac.kr) : 2019/02/06, 2020/04/30
#
#
#         DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.
#
#         INPUTS:
#         1. Plink dataset (*.bed/bim/fam)
#         2. Reference dataset (*.bgl.phased.vcf.gz(Beagle 4.1))
#
#         DEPENDENCIES: (download and place them in the 'dependency/' folder.)
#         1. PLINK (v1.9)  (Will not work with older Plink 1.07)
#         2. Beagle (v4.1) (Need to rename java executable as 'beagle.jar')
#         3. vcf2gprobs.jar (Beagle utility for generating a Beagle v3 genotypes probability file from a Beagle 4.1 vcf file with GT field data)
#         4. [Optional] If genetic_map_file argument is specified, PLINK format genetic map on cM scale
#             (plink.chr6.GRCh36.map, downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)
#
#
#     #################################################################################################
#                                      '''))
#
#
#     parser._optionals.title = "OPTIONS"
#
#     parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')
#
#     parser.add_argument("--target", "-t", help="\nThe prefix of Target PLINK SNP data to impute.(.bed/.bim/.fam)\n\n", required=True)
#     parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)
#     parser.add_argument("--reference", "-ref", help="\nThe prefix of reference panel for imputation.\n\n", required=True)
#
#     parser.add_argument("--tolerated-diff", help="\nTolerated diff (default : 0.15).\n\n", default=0.15)
#     parser.add_argument("--dependency", help="\nSpecify a folder for external software.\n\n", default='dependency/')
#
#     # Beagle(v4).
#     parser.add_argument("--mem", help="\nJava Memory requried for Bealge4.1. (ex. 2g)\n\n", default="2g")
#     parser.add_argument("--nthreads", help="\nThe number of threads to use in Bealge4.1. (ex. 2)\n\n", default=1, type=int)
#     parser.add_argument("--niterations", help="\n'niterations' parameter in Beagle4.1 (default: 5).\n\n", default=5)
#
#
#     ##### <for Test> #####
#
#
#     ##### <for Publication> #####
#
#     args = parser.parse_args()
#     print(args)
#
#     SNP2HLA(args.target, args.reference, args.out, _mem=args.mem, _tolerated_diff=args.tolerated_diff,
#             _p_dependency=args.dependency, _b_nthreads=args.nthreads, _b_niterations=args.niterations)
