# *Isopogon* ddRAD analyses  
Author: B.M. Anderson  

These notes outline the analyses run on ddRAD data to generate results for the paper "[Title]"  
All scripts are assumed to be in an accessible folder for relative calls, indicated here as a directory in home called `~/scripts/`  


# Data
ddRAD enzymes: PstI and MseI  
Illumina NovaSeq 6000, 150 bp single end reads, four indices: ACAGTG, CTTGTA, GCCAAT, GTGAAA (2 files per index)  

Barcodes are per index, so there will be four sets to demultiplex (and four barcode files needed)  
Create a tab-delimited text file (`barcodes.tab`) with columns corresponding to (no header): sampleID seqID index barcode
```s
for index in ACAGTG CTTGTA GCCAAT GTGAAA
do
# grab lines with the index
grep "$index" barcodes.tab > temp
# join the first two columns with a "-" and add a column with the barcode
paste <(cut -f 1,2 temp | tr "\t" "-") <(cut -f 4 temp) > barcodes_"$index".tab
done
rm temp
```
This will produce four barcodes files, one per index  


# Assembly in ipyrad
ipyrad v. 0.9.81 running in a Singularity container (ipyrad installed with miniconda)  
Depending on the high performance computing setup, the method of executing the run will change  
The run was executed on the Zeus cluster on Pawsey (no longer running), using an sbatch script to call the container to execute a Python script that specified the run parameters (`ipyrad_full.py`)  

An initial run showed that 21 samples had too few reads to be retained:
BUX1-06-344560
DAR1-07-344576
FIT1-03-344581
FIT3-05-344594
FIT4-05-344677
OBO1-03-344597
OBO2-05-344685
OBO3-01-344687
OBO3-02-344688
POL1-05-344607
POL1-06-344608
POL2-06-344615
POL6-03-344705
RAV1-03-344626
RAV1-04-344711
RAV1-05-344712
RAV2-01-344713
RAV3-03-344720
ZIA1-02-344640
ZIF1-01-344736
ZIF1-03-344738

Based on preliminary optimisation runs on subsets of the data, the optimal clustering threshold was determined to be 0.92  

The full run therefore used the optimal clustering threshold and was subset to exclude the above 21 samples (excluded from a `samples_full.txt` file containing all samples to run)  
The run used 28 CPUs and took roughly 10 hours  
The command to run was roughly as follows, with the read and barcode files in the working directory
```s
singularity exec -B "$(pwd)":/home/"$USER" ipyrad.sif python ipyrad_full.py 28 samples_full.txt
```
The primary output for downstream analyses is the VCF file: `full.vcf`  
The VCF file is largely unfiltered, keeping all loci found in at least three samples  
Another output, `full.loci`, can be used to access all sites (not just SNPs)  


# Filtering
Datasets were filtered for each analysis differently, and were informed by an initial estimate of error rate based on technical replicates  
During library prep and sequencing, five samples failed and had their spots taken by duplicates of existing samples:
```
Failure	Replacement (duplicating the original)
OBO1-02-344596	FIT2-04-344588
OBO2-03-344602	FIT3-04-344593
POL1-07-344609	POP1-01-344619
POL2-04-344613	BUX3-04-344656
POL5-01-344700	POP2-01-344706 (note: already replicated, so now a triplicate)
```
As a result, those sample labels in the output VCF represent additional duplicates  

After estimating error **and clonality (see below)**, one of each of the technical replicates can be removed (chosen manually based on more recovered loci; the above extra duplicates were all removed; two were removed from the triplicate); these sampleIDs can be put in a `reps.to_remove` text file, one per line for later filtering steps  


## Tiger genotyping error estimate
Create a subset of the VCF for the replicate pairs (one of the reps in each pair is named with a "_R" subscript for sampleID)  
First, create a list of the samples
```s
grep "_R" samples_full.txt | cut -f 1 -d "_" > temp
grep -f temp samples_full.txt | sort > reps.txt
rm temp
```
Manually add the above five failures and additional replicates to `reps.txt` (but not POP2-01-344706 which is already in the file)  

Now filter the VCF for just those samples, keeping SNPs present in all of them, with a minimum mean depth of 10 and only biallelic SNPs using the script `filter_vcf.py`  
Then use the script `single_snp.py` to select a single SNP per locus (max depth or random in the case of ties in coverage)  
```s
grep -v -f reps.txt samples_full.txt > samples.to_remove 
python ~/scripts/filter_vcf.py -o error --mincov 1 --minmd 10 --bial yes -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes error.vcf
mv mod_error.vcf error.vcf
```

Using the technical replicates and the software Tiger (https://bitbucket.org/wegmannlab/tiger/wiki/Home), estimate error rates  
Create an input samples file for Tiger (needs reps assigned to groups)  
```s
mkdir tiger && cd tiger
num=$(($(wc -l < ../reps.txt) / 2))
for ((i=1; i<=num; i++)); do echo $i >> temp && echo $i >> temp; done
echo -e "Sample\tGroup" > rep_groups.txt
paste ../reps.txt temp >> rep_groups.txt
rm temp
```
Manually assign POL5-01-344700 to the POP2-01 group (triplicate)

Now run error estimation and visualise it with the script `Tiger_error_plot.py`  
```s
tiger task=estimateIndReps outname=error vcf=../error.vcf groups=rep_groups.txt
python ~/scripts/Tiger_error_plot.py -o error -m 500 error_errorRates.txt
```
This will create an `error.pdf` file  
Based on the error rate estimates, mean read depth cutoffs of min 10 and max 250 were chosen  


## Clonality
To assess whether samples are indistinguishable from clones, technical replicates can again be used to establish what level of similarity is expected for two identical samples  
This can be based on the proportion of sites called in both that differ  

First, filter the VCF for only ingroup samples and for SNPs found in at least 90% of the samples (to make dataset size more manageable)  
In our naming convention, outgroup samples have a "Z" in their sampleIDs  
The script to assess similarity is `vcf_similarity.py`  
```s
# back in the main directory (not in the tiger folder)
grep "Z" samples_full.txt > samples.to_remove
python ~/scripts/filter_vcf.py -o clones --mincov 0.9 -s samples.to_remove full.vcf
python ~/scripts/vcf_similarity.py -v clones.vcf -o clones
```

The output `clones_comps.txt` can be imported into a spreadsheet program and sorted to manually check  
The output `clones_hist.png` can be useful to see the distribution of comparisons and whether there is a clear break between the replicate comparisions and other comparisons in the dataset  

No other comparisons had difference less than twice the average difference between replicates (by an order of magnitude), so no samples were considered clonal  


## Set 1: ingroup (strict filtering)
For analyses sensitive to missing data (PCA, Structure), the VCF was filtered to only ingroup samples and biallelic SNPs present in at least 90% of samples (one per locus), with a minimum minor allele count of 3  
Create a list of samples that should be excluded (including replicates and outgroups)  
```s
cat reps.to_remove > samples.to_remove
grep "Z" samples_full.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o set1 --mincov 0.9 --minmd 10 --maxmd 250 --mac 3 --bial yes -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes set1.vcf
mv mod_set1.vcf set1.vcf
```


### Set 1 subset 1
To zoom into portions of the PCA and to run Structure analyses within each cluster of samples, SNPs were filtered for just those samples in a similar way to Set 1  
```s
cat reps.to_remove > samples.to_remove
grep "Z" samples_full.txt >> samples.to_remove
grep "RAV" samples.txt >> samples.to_remove
grep "POP" samples.txt >> samples.to_remove
grep "POL" samples.txt >> samples.to_remove
grep "DAR" samples.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o set1_1 --mincov 0.9 --minmd 10 --maxmd 250 --mac 3 --bial yes -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes set1_1.vcf
mv mod_set1_1.vcf set1_1.vcf
```


### Set 1 subset 2
Likewise for a second subset of samples  
```s
cat reps.to_remove > samples.to_remove
grep "Z" samples_full.txt >> samples.to_remove
grep "BUX" samples.txt >> samples.to_remove
grep "OBO" samples.txt >> samples.to_remove
grep "FIT" samples.txt >> samples.to_remove
grep "SPA" samples.txt >> samples.to_remove
grep "CAN" samples.txt >> samples.to_remove
grep "DAR" samples.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o set1_2 --mincov 0.9 --minmd 10 --maxmd 250 --mac 3 --bial yes -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes set1_2.vcf
mv mod_set1_2.vcf set1_2.vcf
```


## Set 2: ingroup (lax filtering)
For analyses less sensitive to missing data (distance network), the VCF was filtered in a similar way to Set 1 but only requiring SNPs be present in 25% of samples (still one per locus), with no restrictions on allele count  
```s
cat reps.to_remove > samples.to_remove
grep "Z" samples_full.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o set2 --mincov 0.25 --minmd 10 --maxmd 250 --mac 1 -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes set2.vcf
mv mod_set2.vcf set2.vcf
```


## Set 3: phylo (concatenation)
To generate alignments, loci were selected after removing one of the outgroups (with too few samples), keeping loci present in at least 50% of samples  
Create a list of samples to remove  
```s
cat reps.to_remove > samples.to_remove
grep "ZIF" samples_full.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o phylo --mincov 0.5 --minmd 10 --maxmd 250 --mac 1 -s samples.to_remove full.vcf
```

From the filtered VCF, loci names can be extracted and used to extract the full alignments from `full.loci` using the script `loci_extract.py` and at the same time removing undesired samples from those loci  
```s
grep -v "#" phylo.vcf | cut -f 1 | uniq | sed 's/RAD_//g' > loci.txt
mkdir loci_extract && cd loci_extract
python ~/scripts/loci_extract.py -l ../loci.txt -s ../samples.to_remove ../full.loci
```

All the loci can be combined into a single alignment using the script `combine_alignments.py`  
In addition, the alignment needs to then be filtered to remove any positions with more than 50% missing data (using the script `clean_alignment.py`)
```s
python ~/scripts/combine_alignments.py -f single *.fasta
python ~/scripts/clean_alignment.py -p 50 combine_out.fasta
mv combine_out_clean.fasta ../concat_loci.fasta
cd .. && rm -r loci_extract
```
The resulting file is `concat_loci.fasta`  


## Set 4: phylo (coalescent)
For running SVDquartets, extract a single SNP per locus from the same loci as were used for concatenation (present in at least 50% of samples)
```s
python ~/scripts/single_snp.py -r yes phylo.vcf
```
This will create the file `mod_phylo.vcf`  


## Set 5: popgen
For population genetic statistics, we will filter out the outgroups plus any populations with fewer than three individuals sampled  
Also filter even more stringently for biallelic SNPs present in at least 95% of samples  
```s
cat reps.to_remove > samples.to_remove
grep "Z" samples_full.txt >> samples.to_remove
grep "CAN2" samples_full.txt >> samples.to_remove
grep "POL5" samples_full.txt >> samples.to_remove
grep "POL6" samples_full.txt >> samples.to_remove
grep "RAV1" samples_full.txt >> samples.to_remove
sort samples.to_remove | uniq > temp
mv temp samples.to_remove
```

Filter
```s
python ~/scripts/filter_vcf.py -o popgen --mincov 0.95 --minmd 10 --maxmd 250 --mac 3 --bial yes -s samples.to_remove full.vcf
python ~/scripts/single_snp.py -r yes popgen.vcf
mv mod_popgen.vcf popgen.vcf
```


# Analyses
Scripts to assess each of the input VCF files for missing data and read depth are `missing_plot.R` and `vcf_depth.py`  

Population number is based on the west-east ordering of the populations
```
1	DAR3
2	CAN1
3	DAR2
4	DAR1
5	CAN2
6	SPA2
7	SPA1
8	BUX3
9	BUX1
10	BUX2
11	SPA3
12	SPA4
13	SPA5
14	OBO1
15	OBO2
16	FIT4
17	OBO3
18	FIT2
19	FIT1
20	POP2
21	POL5
22	RAV1
23	RAV2
24	RAV3
25	POL6
26	FIT3
27	POL4
28	POP1
29	POL3
30	FIT5
31	POL1
32	POL2
```


## PCA
With datasets 1, 1a and 1b as input, the PCA can be run interactively using `pca.rmd`  
Copy or link the `set1.vcf`, `set1_1.vcf` and `set1_2.vcf` into the working directory where the PCA is run  


## Distances
With dataset 2 as input, run the script `distances.R` to generate a distances matrix  
Use the option `-d M` for MATCHSTATES distances from the package `pofadinr`  
```s
Rscript ~/scripts/distances.R -d M -o set2 -v set2.vcf
```
This will produce a number of files, but the key output is `set2_distMATCHSTATES.nex`, which can be input into SplitsTree4 v. 4.17.1 to create a NeighborNet network  


## Structure
Using dataset 1 subsets 1a and 1b, convert the VCFs to the input format appropriate for Structure, including the generation of population mapping files to order groups (based on results of PCA and Distances, and using geographic location)  

Set 1 subset a
```s
grep "#CHROM" set1_1.vcf | cut -f 1-9 --complement | tr -s "\t" "\n" > temp
paste temp <(cut -f 1 -d "-" temp) > str_pops1a.tab && rm temp
sed -i 's/BUX[1-3]$/1/g' str_pops1a.tab
sed -i 's/CAN[1-2]$/2/g' str_pops1a.tab
sed -i 's/SPA[1-5]$/2/g' str_pops1a.tab
sed -i 's/OBO[1-3]$/3/g' str_pops1a.tab
sed -i 's/FIT[1-2]$/4/g' str_pops1a.tab
sed -i 's/FIT4$/4/g' str_pops1a.tab
sed -i 's/FIT3$/5/g' str_pops1a.tab
sed -i 's/FIT5$/5/g' str_pops1a.tab
```

Set 1 subset b
```s
grep "#CHROM" set1_2.vcf | cut -f 1-9 --complement | tr -s "\t" "\n" > temp
paste temp <(cut -f 1 -d "-" temp) > str_pops1b.tab && rm temp
sed -i 's/POL[1-6]$/1/g' str_pops1b.tab
sed -i 's/POP2$/1/g' str_pops1b.tab
sed -i 's/POP1$/2/g' str_pops1b.tab
sed -i 's/RAV[1-3]$/3/g' str_pops1b.tab
```

Use the mapping files and generate the Structure input files with the script `vcf_to_structure.py`  
```s
python ~/scripts/vcf_to_structure.py -p str_pops1a.tab set1_1.vcf
mv set1_1.vcf.str set1_1.str
python ~/scripts/vcf_to_structure.py -p str_pops1b.tab set1_2.vcf
mv set1_2.vcf.str set1_2.str
```
`Converted a VCF file with 83 samples and 9963 SNP loci to Structure format with 83 samples and 9963 SNP loci`  
`Converted a VCF file with 43 samples and 7191 SNP loci to Structure format with 43 samples and 7191 SNP loci`  

Use the `mainparams` and `extraparams` files from Structure (see https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html), and alter specific lines to (as appropriate for the input files just generated):
```s
# mainparams
BURNIN	100000
NUMREPS	100000
NUMINDS	83
NUMLOCI	9963
MARKERNAMES	0

# for 1b
NUMINDS	43
NUMLOCI	7191

#extraparams
UPDATEFREQ	1000
RANDOMIZE	0
```

Running Structure can be done on a cluster, submitting multiple jobs, one per K, each run using the input file `set1_1.str` or `set1_2.str` and the corresponding two params files  
In this case, Structure was run for 40 replicates (in parallel) per K from K1 to K7 (14 20-core jobs) for each of the two sets  
For example, using a batch script:  
```s
for ((kval=1; kval<=7; kval++)); do
sbatch ../structure.sbatch -k $kval -s set1_1.str &&
sbatch ../structure.sbatch -k $kval -s set1_1.str -r 21
done
```
Each rep in each job called Structure like so:
```s
structure -m mainparams -e extraparams -K "$k_val" -i "$str_file" -o output_"$k_val"_"$rep" -D "$rep""$k_val""$RANDOM"
```

The runs took from roughly 1 to 6 hours each on Pawsey (Zeus)  

For reducing file size, remove all the extraneous "Locus" information from the output files  
(from within the output folder, having files named in the pattern "output_x_x_f")
```s
for file in *_f; do
sed -n '/Locus/q;p' $file > temp && mv temp $file
done
```

Evaluate the likelihoods of the replicates and keep the top 10 per K
```s
# grab likelihoods from the output files
for K in {1..7}; do
for rep in {1..40}; do
grep "Estimated Ln Prob" output_"$K"_"$rep"_f | cut -f 2 -d "=" | sed "s/^[[:space:]]*/$K\t/" >> k${K}_likes.txt
done
done

for K in {1..7}; do
# number the likelihoods to keep track of which reps they refer to
num=$(wc -l < k"$K"_likes.txt)
for ((i=1; i<=num; i++)); do
echo $i >> temp.txt
done
paste k"$K"_likes.txt temp.txt > temp_lines.txt && rm temp.txt
# sort the likelihoods and keep the top 10
sort -k2 -nr temp_lines.txt | cut -f 3 | head -n 10 > temp_select.txt
# grab those reps and lines from the likelihood files for plotting
for row in $(cat temp_select.txt); do
cp output_"$K"_"$row"_f select_"$K"_"$row"_f
done
cut -f 1,2 <(sort -k2 -nr temp_lines.txt) | head -n 10 > select_k"$K"_likes.txt
done
rm temp_lines.txt temp_select.txt
```

Evaluate the best K based on the Evanno method using the script `bestK_Evanno.py`  
```s
cat select*likes.txt > all_likes_selected.txt
python ~/scripts/bestK_Evanno.py -l all_likes_selected.txt -o set1_1_bestK_selected
```

Upload the top 10 per K to CLUMPAK (http://clumpak.tau.ac.il/index.html) to run with default settings
```s
zip results.zip select_*_f
```
Upload `results.zip`  

Download, then create qfiles from the major modes for plotting  
Note: in the download, there should be a summary logfile indicating how many minor modes there are per K  
Wherever the downloaded folder is (has subfolders for each K), store that path in a variable e.g. `output=path/to/output`  
```s
for K in {1..7}; do
cut -f 2 -d ':' ${output}/K\=${K}/MajorCluster/CLUMPP.files/ClumppIndFile.output | awk '{$1=$1; print}' > qfile"$K".txt
done
```

If there are minor modes of interest, they can also be kept as qfiles for plotting  
```s
for K in {1..7}; do
for minor_num in {1..4}; do
if [ -d ${output}/K\=${K}/MinorCluster${minor_num} ]; then
cut -f 2 -d ':' ${output}/K\=${K}/MinorCluster${minor_num}/CLUMPP.files/ClumppIndFile.output | awk '{$1=$1; print}' > qfile"$K"_minor"$minor_num".txt
fi
done
done
```

Plot modes of interest using the script `structure_barplots.py`  
Colours for the bars are specified in a `colours.txt` file, one colour (in hex code) per line  
Note: change the qfile and `-o` arguments for the minor mode files, if desired  
Note: to sort the barplots to a specific sampling order, create a text file with samples (one per line) in the order desired and specify it in the call with `-s`  
```s
for num in {2..5}; do
python ~/scripts/structure_barplots.py -o K"$num" -p ../str_pops1a.tab -q qfile"$num".txt -c colours.txt -s sorting.txt
done
```
Repeat the above for set 1 subset b as well  
The resulting barplots can be combined manually into a single figure in Inkscape  


## ML concatenation
Run the concatenated loci (`concat_loci.fasta`) in IQ-TREE v. 2.1.3, running a search for the best model, 1000 ultrafast bootstrap replicates and Shimodaira-Hasegawa-like approximate likelihood ratio tests  
Here it is (essentially) run on Pawsey with maximum 14 cores (took about 11 hours)  
```s
iqtree -s concat_loci.fasta --prefix phylo --threads-max 14 -T AUTO --seqtype DNA -m MFP --merit BIC --ufboot 1000 --alrt 1000
```

The resulting tree (`phylo.treefile`) can be interactively plotted using `plot_tree.rmd`  
Create an `outgroup.txt` file for the process  
```s
grep ">Z" concat_loci.fasta | cut -f 1 -d " " | sed 's/>//g' > outgroup.txt
```

If wanting to collapse clades, make a `clades.txt` file for clades of samples (search terms) that should be collapsed, one per line, tab separated if more than one search term per clade  

It can sometimes be insightful to look at splits present in the bootstrap trees, though in this case the support is high for most nodes  
Use the `phylo.splits.nex` from the run in SplitsTree4 v. 4.17.1 to create a network, using the weights of the splits (bootstrap occurrences), which can be displayed to show other relationships recovered in the bootstrap trees  


## SVDquartets
Use the script `vcf_to_svd.py` to convert the filtered VCF (`mod_phylo.vcf`) into the correct format (and as a single sequence with ambiguities rather than two haplotypes)  
First, set up the hypothesised species to run a grouping of lineages if desired  
```s
grep "#CHROM" mod_phylo.vcf | cut -f 1-9 --complement | tr -s "\t" "\n" > temp
paste temp <(cut -f 1 -d "-" temp) > spp.tab && rm temp
sed -i 's/BUX[1-3]$/1a/g' spp.tab
sed -i 's/CAN[1-2]$/1b/g' spp.tab
sed -i 's/SPA[1-5]$/1b/g' spp.tab
sed -i 's/OBO[1-3]$/1b/g' spp.tab
sed -i 's/FIT[1-2]$/1b/g' spp.tab
sed -i 's/FIT4$/1b/g' spp.tab
sed -i 's/FIT3$/1c/g' spp.tab
sed -i 's/FIT5$/1c/g' spp.tab
sed -i 's/POL[1-6]$/2a/g' spp.tab
sed -i 's/POP2$/2a/g' spp.tab
sed -i 's/POP1$/2b/g' spp.tab
sed -i 's/RAV[1-3]$/2c/g' spp.tab
sed -i 's/DAR[1-3]$/3/g' spp.tab
sed -i 's/ZIA1$/out/g' spp.tab
```

Use that for creating the input NEXUS file (`mod_phylo.nex`)  
```s
python ~/scripts/vcf_to_svd.py -v mod_phylo.vcf -f seq -s spp.tab
```

Now create batch files (text files called `batch.nex`) for running PAUP (set threads appropriately; here for on a cluster), the first for lineages  
```s
#NEXUS
BEGIN PAUP;
	LOG start file=log.txt;
	SET autoclose=yes;
	exe mod_phylo.nex;
	svdq nthreads=28 evalQuartets=all treeInf=QFM bootstrap=no ambigs=distribute;
	describeTrees;
	saveTrees file=qfm.tre format=Newick;
	svdq nthreads=28 evalQuartets=all treeInf=QFM bootstrap=standard nreps=100 treefile=boots.tre ambigs=distribute;
	describeTrees;
	saveTrees file=bootcon.tre supportValues=nodeLabels format=Newick;
	LOG stop;
	QUIT;
END;
```

and the second for hypothesised species  
```s
#NEXUS
BEGIN PAUP;
	LOG start file=log.txt;
	SET autoclose=yes;
	exe mod_phylo.nex;
	svdq nthreads=28 taxpartition=svdtaxa evalQuartets=all treeInf=QFM bootstrap=no ambigs=distribute;
	describeTrees;
	saveTrees file=qfm.tre format=Newick;
	svdq nthreads=28 taxpartition=svdtaxa evalQuartets=all treeInf=QFM bootstrap=standard nreps=100 treefile=boots.tre ambigs=distribute;
	describeTrees;
	saveTrees file=bootcon.tre supportValues=nodeLabels format=Newick;
	LOG stop;
	QUIT;
END;
```

Each batch file should be moved to its own directory along with a copy of (or link to) the nexus input file  
This is for running PAUP* 4.0a (build 168) for Unix/Linux in a Singularity container, but it gives the rough idea for running PAUP  
```s
singularity exec paup.sif paup batch.nex
```
(this would need to be executed separately for the lineages and species runs, in separate folders)  

The resulting trees can be plotted using the same `plot_tree.rmd`  

For the species run, create an `outgroup.txt` manually with the line `out`  
Plot the `qfm.tre`, then also the `bootcon.tre` to get bootstrap values; if any clades differ, the bootstrap support can be obtained from the log file  

For the lineages run, create an `outgroup.txt` and `clades.txt` for collapsing groups  
```s
cut -f 1 mod_phylo.nex | grep "Z" > outgroup.txt
```
Note that bootstrap support is not in the `qfm.tre` file, so it can be input using the values from the log file or the `bootcon.tre` file (though sometimes the branching will differ from the inferred tree)  
If the grouping is not present in the `bootcon.tre`, use IQTREE to map the bootstrap trees onto the obtained qfm tree like so: 
(note: remove any underscores, as this causes an error; also, for some reason IQTREE reports bootstrap support * 100) 
```s
sed "s/_R/R/g" boots.tre > boots_mod.tre
sed "s/_R/R/g" qfm.tre > qfm_boots.tre
iqtree -sup qfm_boots.tre -t boots_mod.tre
mv boots_mod.tre.suptree qfm_boots.tre
```

If wanting to show site concordance factors in addition to bootstrap values, they can be calculated with IQ-TREE (for the lineages run)  
Note that the output trees need to be adjusted to have branch lengths to run properly (just set to 1 in this case)  
```s
sed "s/,/:1,/g" qfm.tre | sed "s/)/:1)/g" > qfm_mod.tre
iqtree -s mod_phylo.nex -t qfm_mod.tre --scf 100000 -T 8 --prefix concord
```
The output can be extracted into a tree form usable by the `plot_tree.rmd` using the script `concord_to_newick.py`  
```s
python ~/scripts/concord_to_newick.py -t concord.cf.tree.nex -o concord
```
This produces the file `concord_scf.tre` to input as a treefile to `plot_tree.rmd`  


It can be insightful to look at the bootstrap tree splits, which can be calculated using the `boots.tre` files from both runs as input to SplitsTree4 v. 4.17.1 to construct Consensus Networks showing splits present in at least 10 trees and weighted by count 






## Popgen
First, create the sample ID and pop ID text file, tab-separated, one per line  
```s
grep "#CHROM" popgen.vcf | cut -f 1-9 --complement | tr -s "\t" "\n" > temp
paste temp <(cut -f 1 -d "-" temp) > popgen_pops.tab && rm temp
```

For looking at measures of diversity including constant sites, we need to extract the corresponding loci and create an alignment  
(again, removing samples using the same file as was used to filter `popgen.vcf`)
```s
grep -v "#" popgen.vcf | cut -f 1 | uniq | sed 's/RAD_//g' > ploci.txt
mkdir loci_extract && cd loci_extract
python ~/scripts/loci_extract.py -l ../ploci.txt -s ../samples.to_remove ../full.loci
```

All the loci can be combined into a single alignment using the script `combine_alignments.py`  
```s
python ~/scripts/combine_alignments.py -f single *.fasta
mv combine_out.fasta ../popgen_loci.fasta
cd .. && rm -r loci_extract
```
The resulting file is `popgen_loci.fasta`  
Note: change the sample labels in the fasta to match what is expected (remove all but the ID from each entry)
```s
sed -i 's/ concat$//g' popgen_loci.fasta
```

Now, evaluate popgen stats for the VCF and for the fasta files with the script `popgen_stats.R`  
This includes running Fst calculations for the VCF file  
```s
Rscript ~/scripts/popgen_stats.R -o popgen -s popgen_pops.tab -v popgen.vcf -f popgen_loci.fasta -t
```

Pairwise Fst values can be plotted in a heatmap with the script `fst_heatmap.R` and a text file (`temp_pops.txt`) with populations, one per line, in the order desired for plotting (optional)  
```s
Rscript ~/scripts/fst_heatmap.R -o popgen -f popgen_StAMPP_Fst.txt -p temp_pops.txt -c mint
```
This will create a pdf and a png called `popgen_Fst`  


### Isolation by distance
We can use the output Fst comparisons and population localities to assess signals of isolation by distance within groups of interest with the script `isolation_by_distance.R`, which can run a Mantel test as well as a linear regression  
Create a population localities file (`pop_loc.tab`) having a single population, lat, and long separated by tabs, one per line  
Also select the populations of interest in a `temp_pops.txt` file (one per line), then run the script like so:  
```s
Rscript ~/scripts/isolation_by_distance.R -o set1 -d popgen_StAMPP_Fst.txt -l pop_loc.tab -m -r -s temp_pops.txt
```

Repeat this for each of the various sampling combinations (new `temp_pops.txt` files)  
Mantel test statistics and significance will be output to the screen for each run  
A pdf of the linear regression plots will be created with the output prefix specified (e.g. `set1_IBD.pdf`)
