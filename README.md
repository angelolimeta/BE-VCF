# BE-VCF
Variant calling for base editors in yeast

**Brief overview of the key analysis steps**
 
- Filter and trim low quality base calls from the raw .fastq data (fastp).
- Index the reference fasta sequence (Samtools) and align the quality filtered reads to it (burrows-wheeler aligner). This results in a .bam file containing the reads aligned to the reference sequence.
- Pileup the .bam file (bcftools). Basically, at each position of the reference sequence, we compile a list of sequences in the reads that differ from the reference sequence. These piled-up sequences represent SNPs or indels.
- Parse the pileup data into the tabular vcf format (bcftools). This is simply done so that the data becomes human-readable and easy to work with in R. The vcf file is simply a tabular data file containing information of SNPs, Indels and sequencing depth at each position in the reference sequence. It contains this information on a per-sample basis.
- Calculate the normalized mutation frequency at each base. So, here, for each position in the reference sequence we count the number of non-reference reads and divide this by the total amount of reads mapped to this particular position. This is done to correct for differing sequencing depth across samples.
- Calculate the mutation enrichment at each base. This is done by simply dividing the mutation normalized frequency in a sample by its corresponding T0 sample. For example, lets say sample T0_dCasAID_gr1 and T1_dCasAID_gr1 has a fraction of 0.01 and 0.02 non-reference bases at position 100, respectively. The mutation enrichment for position 100 in the T1 sample is therefore 0.02 / 0.01 = 2 at that particular position.
- Plot mutation enrichment for all positions in the reference sequence across all samples and compare.
