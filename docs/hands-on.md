# Merqury Hands-On

Get a "Large" instance: 8 cores, 16 G RAM, 50 GB Storage

## Learning objectives

- Count kmers
- Basic histogram and numbers in the kmer db
- Run Merqury on various assemblies
- Understand Merqury outputs: QV, spectra-cn, spectra-asm, phase block, switches

## Where is the data?

- Download https://genomeark.s3.amazonaws.com/trainingmaterials/BGA2023/Merqury/data.tar.gz under /workspace/ and unzip, untar.
- We assume all the input is under `/workspace/data/`.
- Parental mers: generated from the parental strain assemblies of the falcon-unzip paper ([Chin et al., Nat. Methods, 2016]([10.1038/nmeth.4035](https://doi.org/10.1038%2Fnmeth.4035))).
- Assembly: generated with verkko, using 50x downsampled HiFi and ONT R10 reads from a synthetic Col-0 and Cvi-0 of the A. thaliana centromere paper ([Wlodzimierz et al., Nature, 2023]([https://doi.org/10.1038/s41586-023-06062-z](https://doi.org/10.1038/s41586-023-06062-z))) with and without the full Col-0 and Cvi-0 data as parental genomes.

## Counting k-mers

Let’s count the kmers in one assembly. Counting operation works in the same way for reads (`.fq.gz`). For the sake of time, we are just doing it for one assembly.

```bash
mkdir -p test && cd test
ln -s ../data/asm/triocanu_clr/f1_triocanu_col.fasta
meryl count k=21 f1_triocanu_col.fasta output test.k21.meryl
```

## Inspect child’s read histogram

Now, let’s inspect the child’s read kmer set.

```bash
cd ../data/f1/
meryl statistics F1.k21.meryl | head
```

<details><summary>💡 How many distinct kmers are in the child’s read?</summary>

```bash
meryl statistics F1.k21.meryl | head

Found 1 command tree.
Number of 21-mers that are:
  unique              220768892  (exactly one instance of the kmer is in the input)
  distinct            **369212673**  (non-redundant kmer sequences in the input)
  present            4426650250  (...)
  missing         4397677298431  (non-redundant kmer sequences not in the input)

             number of   cumulative   cumulative     presence
              distinct     fraction     fraction   in dataset
frequency        kmers     distinct        total       (1e-6)
--------- ------------ ------------ ------------ ------------
```

  Why is it much larger than the genome?
</details>

Generate the histogram:

```bash
meryl histogram F1.k21.meryl > F1.k21.hist
```

Download `F1.k21.hist`, and let [GenomeScope](http://qb.cshl.edu/genomescope/genomescope2.0/) tell you about this genome.

<details><summary> 💡 What’s the genome size? Level of heterozygosity? </summary>

Here is a pre-run GenomeScope result: [http://genomescope.org/genomescope2.0/analysis.php?code=fPC4nXMZtsJF5o5zCPCZ](http://genomescope.org/genomescope2.0/analysis.php?code=fPC4nXMZtsJF5o5zCPCZ)

![Untitled](https://github.com/BGAcademy23/merqury/assets/12814549/ca58d6a1-7128-4416-98b4-88403b017a15)

Answer: Haploid genome size ~153 Mb. Heterozygisoty level 0.991%.
</details>

## Merqury with no parental data

Make a new directory, and link the data needed.

```bash
mkdir -p /workspace/hands-on/verkko_asm/ && cd /workspace/hands-on/verkko_asm/

# Assembly
ln -s /workspace/data/asm/verkko_hifi/f1_verkko_hifi.fasta

# Read kmer db
ln -s /workspace/data/f1/F1.k21.meryl

# See how to run merqury
$MERQURY/merqury.sh
Usage: merqury.sh [-c] <read-db.meryl> [<mat.meryl> <pat.meryl>] <asm1.fasta> [asm2.fasta] <out>
	-c		: [OPTIONAL] input meryl databases are homopolymer compressed, while asm.fasta isn't
			  This option is not supported for trio-based analysis. Use pre-compressed assemblies with no -c option.
	<read-db.meryl>	: k-mer counts of the read set
	<mat.meryl>	: k-mer counts of the maternal haplotype (ex. mat.hapmer.meryl)
	<pat.meryl>	: k-mer counts of the paternal haplotype (ex. pat.hapmer.meryl)
	<asm1.fasta>	: Assembly fasta file (ex. pri.fasta, hap1.fasta or maternal.fasta)
	[asm2.fasta]	: Additional fasta file (ex. alt.fasta, hap2.fasta or paternal.fasta)
	*asm1.meryl and asm2.meryl will be generated. Avoid using the same names as the hap-mer dbs
	<out>		: Output prefix
Arang Rhie, 2022-09-07. arrhie@gmail.com
```

### Run Merqury (5-10 min)

```bash
$MERQURY/merqury.sh F1.k21.meryl/ f1_verkko_hifi.fasta verkko_asm
```

<details><summary> ⏳ While waiting… What kmer size should I use? </summary>

Merqury has a simple script that computes the best k-size based on [Fofanov et al., Bioinformatics (2004)]([10.1093/bioinformatics/bth266](https://doi.org/10.1093/bioinformatics/bth266)).

```
$MERQURY/best_k.sh
Usage: ./best_k.sh [-c] <genome_size> [tolerable_collision_rate]
  -c         : [OPTIONAL] evaluation will be on homopolymer compressed genome. EXPERIMENTAL
  genome_size: Haploid genome size or diploid genome size, depending on what we evaluate. In bp.
  tolerable_collision_rate: [OPTIONAL] Error rate in the read set. DEFAULT=0.001 for illumina WGS
See Fofanov et al. Bioinformatics, 2004 for more details.
```

Depending on the error rate in the reads, the k could be smaller or larger. The only thing to remember is that smaller kmer have a higher collision rate, meaning there is a higher chance the repeats will be less distinguishable. I recommended using k=21 for reasonable sized genomes (~5 Gb) when evaluating CLR-based assemblies. I’d like to note that the k given here is the lower bound of the k-size. Depending on the assembly quality and the bp accuracy of the sequencing read set, it is recommended to increase the k size, e.g. k=31.

For a. thaliana, the haploid genome size was ~153 Mb - so the best k for evaluating a haploid assembly is

```
$MERQURY/best_k.sh 153000000
genome: 153000000
tolerable collision rate: 0.001
18.5766
```

For diploid assembly, simply double the genome size as a rough estimate:

```
$MERQURY/best_k.sh 306000000
genome: 306000000
tolerable collision rate: 0.001
19.0766
```

We are using k=21 throughout this exercise, which is higher than the suggested 19. The higher the k-mer size, the more conservative the QV becomes.
</details>

## What’s in the output?

### QV

```bash
cat verkko_asm.qv
# f1_verkko_hifi	620933	304044764	40.1169	9.73443e-05
```
<details><summary> 💡 What is the QV? </summary>

QV is `40.1`. There are `620933` error kmers, kmers never seen in the reads. The assembly has `304044764` kmers. This matches roughly the diploid genome size.

</details>

### Spectra-cn plots

See `verkko_asm.f1_verkko_hifi.spectra-cn.ln.png`. This is the ‘unstacked’, line plot of the kmer spectrum.

<details><summary>  💡 verkko_asm.f1_verkko_hifi.spectra-cn.ln.png </summary>
  
![verkko_asm f1_verkko_hifi spectra-cn ln](https://github.com/BGAcademy23/merqury/assets/12814549/6255fbc3-3fbe-4e47-aac5-3caab1feea8a)

</details>

<details><summary> 💡 How does the assembly spectrum look like? </summary>

The kmers seen once (red 1) in the assembly overlap the haploid, 1-cp peak. Likewise, kmers seen twice (blue 2) in the assembly overlap the 2-cp peak. The 3-cp and 4-cp peak are also observed in the expected range, albeit it is hard to clearly distinguish the peaks. In overall, the assembly seems to contain copy numbers as expected, within the 1-4 copy range.

The number of kmers only seen in the assembly (error kmers) are shown as a bar at 0. These are the kmers never seen in the read, thus likely representing a consensus error.

The read-only curve shows kmers present in the read, but never found in the assembly. These kmers are usually occur at low-frequency across the read set, and are likely errors in the sequenced read. One should be alerted when the black line has a ‘bump’ at 1-cp peak or higher, as those are indicating missing sequences in the assembly. Completeness is measured using these kmers.

</details>

### Spectra-asm plots

This time, we are interested in the assembly spectrum, regardless of the kmer counts in the assembly. This helps visualize the missing portion of the assembly.

<details><summary> 💡 See verkko_asm.spectra-asm.ln.png </summary>

![verkko_asm spectra-asm ln](https://github.com/BGAcademy23/merqury/assets/12814549/64e25393-0ad5-44a8-a1c4-b836d36c2510)

</details>

### Completeness

Merqury uses heuristics to define ‘reliable’ kmers. Ignoring the low-frequency kmers, this time, Merqury computes how many are observed in the assembly out of the total reliable kmers.

```bash
cat verkko_asm.completeness.stats 
f1_verkko_hifi	all	129890454	131495725	98.7792
```

<details><summary> 💡 What’s the threshold used to compute completeness? </summary>

See `logs/verkko_asm.spectra-cn.log`.

```bash
...
Filter out kmers <= 5
...
```

This means we are using kmers seen more than 5 times in the reads as ‘reliable’. Note there is always a chance a ‘read-only’ kmer to be seen in the ‘reliable’ set, regardless of the sequencing depth, due to sequencing errors and biases, which could probabilistically occur at a higher frequency.
</details>

## Merqury with parental data

### Genrate hap-mers

This time, we will utilize the parental kmers. The pre-computed meryl dbs are available under `/workspace/data/parental/`. It is important to run `[hapmers.sh](http://hapmers.sh)` before using hapmers.

```bash
# Make a new directory, 'hapmer'.
cd ../
mkdir -p hapmers

# Link data
ln -s ../../data/f1/F1.k21.meryl           # F1 read kmer
ln -s ../../data/parental/Col-0.k21.meryl/ # Col-0 haplotype genome asm's kmer
ln -s ../../data/parental/Cvi-0.k21.meryl/ # Cvi-0 haplotype genome asm's kmer

# See how to run hapmers.sh
$MERQURY/trio/hapmers.sh
```

<details><summary> 💡 How do I generate my hapmers? </summary>

Because we are using kmers from haplotype genome’s assemblies, and not directly from the reads, use `-no-filt` option. This is an unusual case, and only done because we are limited from space and time for this session. If possible, count kmers from the parental reads directly and run `hapmers.sh` without the `-no-filt` option.

```bash
# Run hapmers.sh
$MERQURY/trio/hapmers.sh Col-0.k21.meryl/ Cvi-0.k21.meryl/ F1.k21.meryl/ -no-filt
```

</details>

<details><summary> 💡 What is the intherited, hapmer spectrum? </summary>

See `inherited_hapmers.ln.png`.

![inherited_hapmers ln](https://github.com/BGAcademy23/merqury/assets/12814549/520328df-cb9e-44c4-8df3-63d5e46ca2fc)

Note the roughly 1:1 portion of the hapmers, which nicely covers most of the 1-cp peak. We see this more often in highly diverged sub-species, with minimal shared genomic content between the parents. Compare this to a hapmer obtained from a human:

![inherited_hapmers ln 1](https://github.com/BGAcademy23/merqury/assets/12814549/653fa2b0-7933-4afc-8400-df3ef48f67ac)

Here, we see a lot less of the mat pat hapmers, and a larger amount of shared kmers in the 1-cp region (10~50x kmer multiplicity). This is likely because the child inherited kmers from the parental genome, which both heterozygous haplotypes are present in the parental genome (Mat had AB and Pat had AB, Child inherited AB) or are not distinguishable (Mat had AA and Pat had AB, Child inherited AB - in this case, A is not distinguishable, only B is.).

</details>

### Running Merqury with hapmers on a haplotype phased, diploid assembly (~20 min)

```bash
cd ../
mkdir -p verkko_dip && cd verkko_dip

ln -s ../../data/f1/F1.k21.meryl/
ln -s ../hapmers/Col-0.k21.hapmer.meryl/
ln -s ../hapmers/Cvi-0.k21.hapmer.meryl/
ln -s ../../data/asm/verkko_trio_hifi/f1_verkko_col.fasta 
ln -s ../../data/asm/verkko_trio_hifi/f1_verkko_cvi.fasta

$MERQURY/merqury.sh F1.k21.meryl/ Col-0.k21.hapmer.meryl/ Cvi-0.k21.hapmer.meryl/ f1_verkko_col.fasta f1_verkko_cvi.fasta verkko_dip

read: F1.k21.meryl

Haplotype dbs provided.
Running Merqury in trio mode...

hap1: Col-0.k21.hapmer.meryl
hap2: Cvi-0.k21.hapmer.meryl
asm1: f1_verkko_col.fasta
asm2: f1_verkko_cvi.fasta
out : verkko_dip

Get spectra-cn plots and QV stats

Get blob plots

Get haplotype specfic spectra-cn plots

Get phase blocks
...
```

### Check QV, completeness, haplotype switches

```bash
cat verkko_dip.qv
f1_verkko_col	277023	148586323	40.513	8.88594e-05
f1_verkko_cvi	213105	154022987	41.8092	6.59288e-05
Both	490128	302609310	41.1246	7.71868e-05
```

This time, we have QV for each and both assemblies combined.

```bash
cat verkko_dip.completeness.stats 
f1_verkko_col	all	108663696	131495725	82.6367
f1_verkko_cvi	all	108589169	131495725	82.58
both	all	129889048	131495725	98.7782
f1_verkko_col	Col-0.k21.hapmer	21749217	21760542	99.948
f1_verkko_col	Cvi-0.k21.hapmer	322137	21623509	1.48975
f1_verkko_cvi	Col-0.k21.hapmer	579708	21760542	2.66403
f1_verkko_cvi	Cvi-0.k21.hapmer	21401749	21623509	98.9744
both	Col-0.k21.hapmer	21751138	21760542	99.9568
both	Cvi-0.k21.hapmer	21610226	21623509	99.9386
```

The first three lines show completeness (%) per assembly and both, as the portion of reliable kmers found in the reads. The rest are shown per-hapmers.

f1_verkko_col	Col-0.k21.hapmer - here we expect only Col-0 kmers to be found. However, the next line, f1_verkko_col	Cvi-0.k21.hapmer indicates there are 1.49% of the Cvi-0 hapmers found in the col assembly, which is likely a switch error. Likewise, f1_verkko_cvi	Col-0.k21.hapmer shows unexpected Col-0 hapmers found in the cvi assembly.

both lines at the bottom shows the haplotype completeness regardless of each assembly.

```bash
cat verkko_dip.*.switches.txt
verkko_dip.f1_verkko_col.100_20000 switch error rate (%) (Num. switches / Total markers found): 261481	25730538	1.01623%
verkko_dip.f1_verkko_cvi.100_20000 switch error rate (%) (Num. switches / Total markers found): 661467	25227589	2.622%
```
Merqury automatically tries to find locally phased blocks given the marker presence. By default, Merqury allows up to 100 switches within 20kb, roughly translated into allowing maximum of 0.5% of short-range switches within a minimum of ~20kb window. Once the phase blocks are defined, the num. observed switches (unexpected hapmers) are counted for calculating the switch error rate (%). Here, we see col assembly has switch error rate of 1.0%, while the cvi assembly has slightly higher, 2.6%. [Find more details here](https://github.com/marbl/merqury/wiki/3.-Phasing-assessment-with-hap-mers#5-phased-block-statistics-and-switch-error-rates).

### Visualizing phase switches

The most intuitive way in Merqury to visualize phase switches is via blob-plot. See `verkko_dip.hapmers.blob.png`.

![verkko_dip hapmers blob](https://github.com/BGAcademy23/merqury/assets/12814549/340a0c8b-6700-4ca6-8d3a-5d6c9c7218fe)

Each circle represents a scaffold, with the circle size relative to the length. The better the assembly is haplotype resolved, we see the circles plotted along the X and Y axis, which indicates the num. of hapmers found in each scaffold.

There are a few off-axis circles, indicating switch error in the assembly.

Next, we can check how large the phased blocks are. See `verkko_dip.f1_verkko_cvi.block.N.png`.

![verkko_dip f1_verkko_cvi block N](https://github.com/BGAcademy23/merqury/assets/12814549/b5bde8be-284a-4724-808c-5aee3a387697)

We see there are a few <~1 Mb blocks consisting of Col hapmers. Where are they in the assembly?

### 💡 Let’s visualize hapmers on IGV

Concatenate each hapmer track to easily visualize on the diploid assembly.

```bash
cat *Col-0.k21.hapmer.wig > verkko_dip.Col-0.k21.hapmer.wig
cat *Cvi-0.k21.hapmer.wig > verkko_dip.Cvi-0.k21.hapmer.wig
```
### On IGV

Under `Genomes`, `Load Genome from URL`: [`https://genomeark.s3.amazonaws.com/trainingmaterials/BGA2023/Merqury/hands-on/verkko_dip.fasta`](https://genomeark.s3.amazonaws.com/trainingmaterials/BGA2023/Merqury/hands-on/verkko_dip.fasta)

Download the verkko_dip.C*-0.k21.hapmer.wig files, or load from URL:

- Under `File`, `Load from URL`:
    - https://genomeark.s3.amazonaws.com/trainingmaterials/BGA2023/Merqury/hands-on/verkko_dip/verkko_dip.Col-0.k21.hapmer.wig
    - https://genomeark.s3.amazonaws.com/trainingmaterials/BGA2023/Merqury/hands-on/verkko_dip/verkko_dip.Cvi-0.k21.hapmer.wig

Also create the error kmer track from *_only.wig:

```bash
cat *_only.wig > verkko_dip.error.wig
```

Again, load it directly or load from URL:

- https://genomeark.s3.amazonaws.com/trainingmaterials/BGA2023/Merqury/hands-on/verkko_dip/verkko_dip.error.wig

💡Tip: Give colors as you’d like for each track.
<img width="1274" alt="image" src="https://github.com/BGAcademy23/merqury/assets/12814549/fe55e446-ffea-4bf8-b089-fa9675ed58a8">

