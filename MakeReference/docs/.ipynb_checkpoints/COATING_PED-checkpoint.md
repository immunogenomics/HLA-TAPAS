# COATING_PED
------

### [ Introduction ]

Currently, HLA allele information is majorly investigated and distributed by IMGT-HLA, and the name of allele is named based on the nomenclature rule which this organzation introduced. This naming system is *de facto* standard for HLA alleles and we needed to generalize our framework to work based on this system.

More specifically, based on the system, an allele name consists of (1) Gene(with '*
(asterisk)' as separator), (2) Fields(1 to maximum 4), (3) optional Suffix. Some alleles will have less than 4-field, but the current standard number of fields for an allele is 4. From now on, i'll use "4-field" as "standard" implicitly.

So far, however, a lot of software and researches have used 2-field allele names including **SNP2HLA** and **MakeReference**. For these software to be compatible with current standard 4-field name system, I introduced this module as to transform previously used 2-field (or 1-field) allele to current standars 4-field allele.  

Furthermore, In IMGT-HLA allele naming system, there are two more ways to classify HLA alleles. These are 'P-group' and 'G-group' classification. It is created to classify alleles based on the major exon information(link). I also get this module to transform given alleles based on these classification rule.

### [ Description ]

Breifly, this module takes (1) .ped file, (2) .iat file, (3) Output format and generates transformed .ped file with requested output format. Here, extension .iat means "I"ntegrated "A"llele "T"able file.  

Grouping된거에서 첫 번째꺼 리턴해주는거라고 설명해줄 것.

### [ Example ]

#### (Case 1): From Standard ped file with various format to 4-field standard names.

I'll assume that users are in "makereference_recode_v2/" directory.

```console

python3 ./src/COATING_PED.py \
        -ped ./data/COATING_PED/DummyPED.standard.100.ped \
        -iat ./data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat \
        -o ./Transformed.standard.100 \
        --4field

```

You don't need to specify '.ped' extension in '-o' option.  

Let's check before and after in case of HLA-A, B, C. Before conducting transformation, the first 5 rows in 'DummyPED.standard.100.ped' are,

A1 | A2 | B1 | B2 | C1 | C2
---|----|----|----|----|----
A*1117|02|B*27|B*35:219|0236|0891
0283N|02:385|38:59|35368|C*02:98|03:02
A*01:138|A*26:42|4212|B*4001|15:109|C*07461
A*30:02|02|B*44261|35:209|C*03|07:258
A*01:01|A*01:79|B*35296|0702|06:55|04:155

After transformation, the equivalent part in the output 'Transformed.standard.100.4field.ped' would be,

A1 | A2 | B1 | B2 | C1 | C2
---|----|----|----|----|----
11:17|02:01:01:01|27:01|35:219|02:36:01|08:91
02:83N|02:385|38:59|35:368|02:98|03:02:01
01:138|26:42|42:12|40:01:01|15:109|07:461
30:02:01:01|02:01:01:01|44:261|35:209|03:02:01|07:258
01:01:01:01|01:79|35:296|07:02:01:01|06:55|04:155

Gene caption("A\*", "B\*", etc.) won't be included in default.

#### (Case 2): From 4-field standard names to P or G-group allele names.  

```console

python3 ./src/COATING_PED.py \
        -ped data/COATING_PED/DummyPED.standard.100.ped \
        -iat data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat \
        -o ./Transformed.standard.100 \
        --G-group

```

This command will generate a new ped file with G-group allele names. Input ped file is same as the one in (Case 1).

A1 | A2 | B1 | B2 | C1 | C2
---|----|----|----|----|----
A*1117|02|B*27|B*35:219|0236|0891
0283N|02:385|38:59|35368|C*02:98|03:02
A*01:138|A*26:42|4212|B*4001|15:109|C*07461
A*30:02|02|B*44261|35:209|C*03|07:258
A*01:01|A*01:79|B*35296|0702|06:55|04:155

Then, those alleles will be transformed into G-group form,

A1 | A2 | B1 | B2 | C1 | C2
---|----|----|----|----|----
11:17|02:01:01G|27:01|35:219|02:36:01|08:91
02:01:01G|02:385|38:59|35:368|02:98|03:02:01G
01:138|26:42|42:12|40:01:01G|15:109|07:461
30:02:01G|02:01:01G|44:261|35:209|03:02:01G|07:258
01:01:01G|01:79|35:296|07:02:01G|06:02:01G|04:155

Likewise, result of transforming 4-field standard names to P-group allele names is,

A1 | A2 | B1 | B2 | C1 | C2
---|----|----|----|----|----
11:17|02:01P|27:01|35:219|02:36P|08:91
0|02:385|38:59|35:368|02:98|03:02P
01:138|26:42|42:12|40:01P|15:109|07:461
30:02P|02:01P|44:261|35:209|03:02P|07:258
01:01P|01:79|35:296|07:02P|06:02P|04:155

#### (Case 3 - 1): From P or G-group alleles to standard 4-field alleles.

Users can also take a ped file with P or G-group alleles as an input to this module. In this case, you need to give a path of the input ped file by using '-ped-Ggroup' or '-ped-Pgroup' instead of '-ped'. For example, 

```console

python3 src/COATING_PED.py \
        -ped-Ggroup data/COATING_PED/DummyPED.Ggroup.Ncap.dc.100.ped \
        -iat data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat \
        -o Transformed.Ggroup.100 \
        --4field

```

This command will transform a ped file with G-group alleles to 4-field standard ped file.  

Before, 

A1 | A2 | B1 | B2 | C1 | C2 
---|----|----|----|----|----
30:02:08|03:01:28|56:48|15:165|05:103:02|06:212
02:01:01G|24:296|58:16:01|15:256|14:21N|07:520
02:200|26:76|07:176|40:341|16:30N|02:12
26:01:37|02:22:01G|07:35|35:66|07:101|04:130
02:674|24:22|15:10:04|15:25:01G|04:184|07:10
32:82|30:11:02|15:01:01G|40:01:01G|03:04:01G|07:29:01

After,  

A1 | A2 | B1 | B2 | C1 | C2 
---|----|----|----|----|----
30:02:08|03:01:28|56:48|15:165|05:103:02|06:212
02:01:01:01|24:296|58:16:01|15:256|14:21N|07:520
02:200|26:76|07:176|40:341|16:30N|02:12
26:01:37|02:22:01:01|07:35|35:66|07:101|04:130
02:674|24:22|15:10:04|15:25:01|04:184|07:10
32:82|30:11:02|15:01:01:01|40:01:01|03:04:01:01|07:29:01

  

Same job taking P-group ped file as input can be conducted.

```console

python3 src/COATING_PED.py \
        -ped-Pgroup data/COATING_PED/DummyPED.Pgroup.Ncap.dc.100.ped \
        -iat data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat \
        -o Transformed.Pgroup.100 \
        --4field

```

Before, 

A1 | A2 | B1 | B2 | C1 | C2 
---|----|----|----|----|----
26:91|02:19|0|35:07|04:224|05:01P
68:34|29:01P|37:01P|27:27|07:532|07:478
02:155|02:307|41:01P|35:04P|06:02P|01:102
02:48|34:09|15:344|15:39P|04:172|08:124
68:113|26:89|35:01P|15:121|05:01P|17:28
03:231P|33:109|27:121|13:23|07:190|07:53

After,  

A1 | A2 | B1 | B2 | C1 | C2 
---|----|----|----|----|----
26:91|02:19|0|35:07|04:224|05:01:01:01
68:34|29:01:01:01|37:01:01:01|27:27|07:532|07:478
02:155|02:307|41:01:01|35:04:01|06:02:01:01|01:102
02:48|34:09|15:344|15:39:01|04:172|08:124
68:113|26:89|35:01:01:01|15:121|05:01:01:01|17:28
03:231:01|33:109|27:121|13:23|07:190|07:53


#### (Case 3 - 2): From P or G-group alleles without double-colon to standard 4-field alleles.

In case that you got a ped file with P or G-group alleles but not without double-colon, where you can't assure digits for those alleles, this module handle this problem. As i said it won't be a prefect guess, however, it will find the most probable and prime candidate.

For example, when you are given a ped file with G-group allele names without double-colon, no digit information to determine fields, then you can get it through,

```console

python src/COATING_PED.py \
       -ped-Ggroup data/COATING_PED/DummyPED.Ggroup.Ncap.Ndc.100.ped \
       -iat data/ClassifyGroups/INTEGRATED_ALLELE_TABLE.imgt3320.iat \
       -o Transformed.Ggroup.100.ndc \
       --4field

```

Before transformation,

A1 | A2 | B1 | B2 | C1 | C2 
---|----|----|----|----|----

After transformation,

A1 | A2 | B1 | B2 | C1 | C2 
---|----|----|----|----|----




(cf) If given ped file has not P or G-group format, then you don't need to consider whether  alleles in ped file have double-colon or not, HLA-gene is captioned or not, etc.  


(cf) If you give an input ped file with P or G-group alleles to the module and order it to make output with same P or G-group format, then it won't do anything because this order is meaningless(Transfomation P-group ped file to P-group or G-group ped file to G-group is meaningless).