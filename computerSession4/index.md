# Computer Session 4: How to calculate polygenic scores
##[Oxford NCRM Summer School](http://www.oxfordsociogenetics.com):
###An introduction to combining social science and molecular genetic research
 

nicola.barban@nuffield.ox.ac.uk
![Ox](/images/ox_brand_special_pos.png)

---

# Outline
1. 
2.
3.

---
# What is a Polygenic score?

* An index that linearly aggregates the estimated effects of individual SNPs on the trait of interest.
*  Can be considered a measure of an individual’s genetic propensity towards a trait.
* In general, a polygenic score for an individual is defined as a weighted sum of a persons genotypes at J loci,

<img src="http://www.sciweavers.org/tex2img.php?eq=PGS_i%3D%5Csum_%7Bj%3D1%7D%5EJ%20x_%7Bij%7Dw_j%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="PGS_i=\sum_{j=1}^J x_{ij}w_j" width="128" height="56" />

 where x is individual *i* 's genotype (0,1,2) at variant *j* and *w* is some sort of weight we can use to construct our scores

---

# What do we need?

* Access to individual-level genotype data of a prediction sample
* GWAS summary statistics from a discovery sample

### For example: [Sociogenome](http://www.sociogenome.com/data/)
![img](/images/sociogenome_data.png)

---


# How polygenic scores differ?

### 1. how to generate the weights *w* ,
### 2. how to determine which  *J* loci to include.

---

# Attention! 

### Remember to use summary statistics from a meta-analysis that **does not include** the genotype data you are using.

[Wray et al. 2013](https://www.nature.com/nrg/journal/v14/n7/pdf/nrg3457.pdf)
![img](/images/pgs_pitfalls_.png)

---


# A very basic polygenic score

* based on on one SNP only **J=1**.
* No weights used, **w=1**.

# rs9930506 SNP, FTO gene.

```bash
grep rs9930506 ../Data/hapmap3-r2_b36_CEU.bim 

16	rs9930506	0	52387966	G	A
 
cat ../Data/FTOscore.txt
rs9930506 G 1

./plink  --bfile ../data/hapmap3-r2_b36_CEU \
		 --score  ../Data/FTOscore.txt sum \
		 --out  rs9930506score
		 	

```

---
# rs9930506 SNP
```bash

cat rs9930506score.log
PLINK v1.90b4.3 64-bit (9 May 2017)
Options in effect:
  --bfile ../Data/hapmap3-r2_b36_CEU
  --out rs9930506score
  --score ../Data/FTOscore.txt sum

Hostname: Nicolas-iMac.local
Working directory: /Users/nicolabarban/Dropbox/sociogenome-scripts/SUMMERSCHOOL2017/computerSession4
Start time: Fri Jun 16 11:17:35 2017

Random number seed: 1497608255
16384 MB RAM detected; reserving 8192 MB for main workspace.
1416121 variants loaded from .bim file.
165 people (80 males, 85 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 112 founders and 53 nonfounders present.
Calculating allele frequencies... done.
Warning: 819 het. haploid genotypes present (see rs9930506score.hh ); many
commands treat these as missing.
Total genotyping rate is 0.996191.
1416121 variants and 165 people pass filters and QC.
Note: No phenotypes present.
--score: 1 valid predictor loaded.
--score: Results written to rs9930506score.profile .

		 	

```


---

# rs9930506 SNP

The option `sum` provides the sum of  Allele count. If not specified, plink gives you the average count.


```bash
 head -20 rs9930506score.profile
 
     FID       IID  PHENO    CNT   CNT2 SCORESUM
    1328   NA06989     -9      2      1        1
    1377   NA11891     -9      2      1        1
    1349   NA11843     -9      2      0        0
    1330   NA12341     -9      2      1        1
    1444   NA12739     -9      2      1        1
    1344   NA10850     -9      2      0        0
    1328   NA06984     -9      2      1        1
    1463   NA12877     -9      2      1        1
    1418   NA12275     -9      2      0        0
   13291   NA06986     -9      2      1        1
    1418   NA12272     -9      2      0        0
    1424   NA10845     -9      2      2        2
    1346   NA10852     -9      2      1        1
   13292   NA07051     -9      2      1        1
    1354   NA12400     -9      2      0        0
   13281   NA12344     -9      2      1        1
    1451   NA12777     -9      2      2        2
    1421   NA12287     -9      2      1        1
    1418   NA10837     -9      2      1        1
 

```

---


#Pruning and Thresholding (P&T). Using PRSice

A more common approach is to calculate a score based on **many** SNPs and use as weights the coefficients from a GWAS discovery.


- PRSice website [link](www.prsice.info)
- PRSice manual [manual](http://www.prsice.info/PRSice_MANUAL_v1.25.pdf)
- Useful paper [link](http://biorxiv.org/content/biorxiv/early/2017/02/05/106062.full.pdf)

---


# An algorithm to calculate polygenic scores

1. Prepare genotype data and reference file (based on summary results) and make them consistent

2. Ensures the markers included in the score are all approximately independent of each other

3.  Omit SNPs whose p-value for association with the phenotype is above a certain threshold.


---


# Prepare files.

* SNPs present in only one data set are removed
* *ambiguous* (A/T or C/G) SNPs are removed

### Ambiguous SNPs

* Since DNA is composed of 2  strands, there is ambiguity over which strand to look at. dbSNP uses the assembled chromosome to establish a plus and a minus strand. S

* If a the genotype is  rs1234(A;A) for a SNP in which the other allele is G, but the reference fike claims that this is a C;T SNP, then logically we can flip the results over and call you a rs1234(T;T). **This is safe and reasonable.**

* Unfortunately if this was instead a SNP where the two alleles are A or T the same flipping logic falls down since both forms  rs1234(A;A) and rs1234(T;T) are possible.

* There are possible solution (mainly based on allele frequency). A reasonable solution when using  **many SNPs** is to **remove ambiguous SNPs.**

----

# Selecting independent SNPs
* Ideally we want to avoid “double-counting” the effects of a causal variants.

* Since variants are in Linkage disequilibrium, the risk is to over count SNPs in more densely genotyped area of the DNA.
* A standard strategy is to use  a **pruning** algorithm that ensures the markers included in the score are all approximately independent of each other

---

# Which SNPs?

### Imputed or genotyped SNPs?
* Depends on the genotyping chip. How many variants have been genotyped?
* Which selection algorithm? Plink uses *clumping* or *pruning*.

### Clumping
* Clumping removes SNPs  in linkage disequilibrium with the local SNP with the smallest GWAS P-value. 

* More info on clumping options can be found [here](https://www.cog-genomics.org/plink/1.9/postproc#clump)

---
#  Starting with PRSice 
Before you can run `PRSice`, the packages that it uses must be installed. 

First, ensure that you have the latest version of `R` downloaded - **if not then update R!**  Next, run the following commands in R:

```R

R
> library(fmsb)
> library(batch)
> library(gtx)
Loading required package: survival
Loading required package: splines
> library(plyr)
> library(ggplot2)

```
If any of these commands generate errors, it means that the corresponding packages have not previously been downloaded and so users must **install these packages** and their dependancies, using the `install.packages()` function in R.

---

# Input data files

PRSice needs the following columns for quantitative traits: ** CHR BP A1 A2   SNP P BETA SE **. 



```bash

gunzip ../Data/EduYears_reference.txt.gz

head ../Data/EduYears_reference.txt

SNP	CHR	BP	A1	A2	EAF	BETA	SE	P
rs11130222	3	49901060	A	T	0.5765	0.026	0.003	4.581e-25
rs2883059	3	49902160	T	C	0.5746	0.026	0.003	7.149e-25
rs3796386	3	49899795	A	G	0.4235	-0.026	0.003	8.883e-25
rs55692411	3	49911155	A	G	0.4254	-0.026	0.002	1.045e-24
rs952594	3	49908023	A	G	0.4254	-0.026	0.002	1.291e-24
rs34654589	3	49911449	C	G	0.5746	0.026	0.002	1.549e-24
rs7613360	3	49916710	T	C	0.3993	-0.026	0.003	2.449e-24
rs11712056	3	49914397	T	C	0.5504	0.025	0.002	7.526e-24
rs148734725	3	49406708	A	G	0.3078	0.027	0.003	1.25e-23

```

--- 


#Prepare genotype data

I use plink to do QC on the genotype file. 
Here I select SNPs with MAF>0.01, missing call rates > 0.05 (--geno), HWE p-value<0.0004
 and individuals with missing genotype> 0.05 (--mind). 

```bash
./plink 	--bfile ../Data/hapmap3-r2_b36_CEU  \
	 						--geno 0.05 \
							--maf 0.01 \
							--mind 0.05 \
							--hwe 1e-4 \
							--make-bed \
							--pheno ../data/EA.phen \
							--out HapMapEA_qc
````
---


### 4. Run PRSice

I specify a new directory for each score. This command calculate scores for the following p-value thresholds: *"5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 0.05,  0.5, 1"*. 

- **base** indicates the summary statistics (prepared)
- **target** indicates the plink file
- **barchart.levels** the p-value thresholds

This is with the so-called hard-called genotype data since the PLINK 1 binary format cannot represent genotype probabilities. PRsice can work also with dosage data. 


``` bash
mkdir EA

R --file=../Software/PRSice_v1.25.R   -q --args  \
  plink  ../plink \
  base ../../Data/EduYears_reference.txt \
  target ../HapMapEA_qc \
  ggfig F \
  score.at.1 T \
  fastscore T \
 barchart.levels    "5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 5e-3, 0.05,  0.5, 1" \
  covary F \
  binary.target F \
  debug.mode T \
	   wd "EA" \
  cleanup T 

```


### 4.b Using LDpred
LDpred needs to be execuded in 3 stages:

1. coord_genotypes.py is used to coordinate the genotype and the summary statistics
2. LDpred.py is used to calculate LD weights (takes a lot of time) and to create score
3. validate.py is used to validate the scores and choose best predicting score based on different fraction of causal variants


This is a test created by charlie to run on the terminal
```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:01:00
#SBATCH --job-name=ldtest
#SBATCH --mail-type=ALL
#SBATCH --mail-user=charlierahal@gmail.com
module purge
module load python/2.7
export PYTHONPATH=$PYTHONPATH:/system/software/linux-x86_64/lib/libplinkio/0.3.1/lib/python2.7/site-packages
python ldpred/coord_genotypes.py --gf=test_data/LDpred_data_p0.001_test_0 --ssf=test_data/LDpred_data_p0.001_ss_0.txt --N=100 --out=bogus.hdf5
python ldpred/LDpred.py --coord=bogus.hdf5 --ld_radius=1 --PS=0.3 --N=100 --out=ldpredoutput
python ldpred/validate.py --vgf=test_data/LDpred_cc_data_p0.001_train_0 --rf=ldpredoutput --out=pgsoutputs
```
This is a file in which I tried to calculate LDpred scores for GIV analysis [file](/ldpred.sh). Still queing in Arcus-b!!!


### 5. Importing scores in Stata and using them for other analysis

You can use something like this:
```Stata
clear all
set more off

 import delimited "../PGS1/PRSice_SCORES_AT_ALL_THRESHOLDS.txt", delimiter(space, collapse) case(preserve) encoding(ISO-8859-1)
 save PGS1_TwinsUK, replace
 
```
### Go to GIV page [link](/GIV.md)


