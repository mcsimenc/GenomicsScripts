# Genomics scripts
Python scripts for manipulating various genomics-related file formats. \
Some definitions are included in multiple scripts so they are as standalone
as possible.

### BED scripts
* `bedIntersect2percentOverlap.py`: output proportion of features participating in an overlap as found by bedtools bedIntersect -wa -wb
* `gff2bed.py`: convert a GFF3 file to a BED format file

### BLAST scripts
* `blast2gff`: convert blastn, blastp, etc. tabular output to GFF3 format
* `blastBestHit.py`: output the highest-scoring hit from blastn, blastp, etc. tabular output (-outfmt 6 or 7)
* `blastFilter.py`: output lines from blastn, blastp, etc. tabular output (-outfmt 6 or 7) where the percent identity satisfies constraints

### Circos scripts
* `coverage2circosLine.py`: calculate average depth of coverage form the output of bedtools genomecov -ibam <bam> -d and output a Circos line track
* `fasta2GCcontentCircosHeatmap.py`: calculate GC content for each window in each sequence in a FASTA file and output a Circos heatmap track
* `fixTrackLabels.py`: replace labels in Circos track file with the integer label from the associated Circos karyotype file
* `gff2circosHeatmap.py`: convert feature coordinates in a GFF3 file to Circos heatmap track format with specified bin size
* `gff2circosTile.py`: convert features in a GFF3 file to Circos tile track format
* `vcfSNPrate2circosLine.py.untested`: takes a VCF file with or without a GFF3 file whose features (genes) coordinates are represented in the VCF file and outputs SNPs rate per gene or a Circos heatmap track of SNP rate/bin size

### FASTA scripts
* `fasta2circosIdeograms`: output sequence lengths as a Circos ideogram file
* `fasta2GCcontentCircosHeatmap.py`: calculate GC content for each window in each sequence in a FASTA file and output a Circos heatmap track
* `fastaExtractSeqs.py`: extract a subset of the sequences in a FASTA file
* `fastaExtractNseqs.py`: extract the first or second or third etc.  n sequences from a FASTA file
* `fastaRenameSeqs.py`: rename FASTA sequence headers according to a mapping of old to new names
* `fastaRenameSeqsByLength.py`: sort FASTA sequences in descending order and rename sequences sequentially
* `fastaSplitSeqs.py`: write a new FASTA file for each sequence in a FASTA file
* `gff2fasta.py`: extract sequences from a FASTA file based on coordinates in a GFF3 file, using the value from a specified key in the GFF3 attributes column

### GFF scripts
* `blast2gff`: convert blastn, blastp, etc. tabular output to GFF3 format
* `gff2bed.py`: convert a GFF3 file to a BED format file
* `gff2circosHeatmap.py`: convert feature coordinates in a GFF3 file to Circos heatmap track format with specified bin size
* `gff2circosTile.py`: convert features in a GFF3 file to Circos tile track format
* `gff2fasta.py`: extract sequences from a FASTA file based on coordinates in a GFF3 file using the value from a specified key in the GFF3 attributes column as the sequence name. Depends on BEDTools and BioPython
* `gff2introns.py`: create a GFF3 with intron features from a GFF3 with gene and exon features or output a list of intron lengths
* `gff3line.py`: contains the GFF3\_class
* `gffAddAttribute.py`: add a key-value pair to the attributes column of a GFF3 file
* `gffFilter.py`: remove or retain GFF3 features on specified scaffolds or with specified values for the ID attribute
* `gffMergeOverlaps.py`: merge overlapping features in a GFF3 file
* `gffRemoveScafPart.py`: remove features in a GFF3 file whose coordinates
* `gffRenameScafs.py`: rename scaffolds in a GFF3 file per a two-column map
* `gffSubset.py`: extracts a subset of a GFF3 file based on values of a chosen attribute key
* `gffv2Exonerate2gff3.py`: convert an Exonerate-generated GFF2 file to GFF3 format

### VCF scripts
* `vcfSNPrate2circosLine.py.untested`: takes a VCF file with or without a GFF3 file whose features (genes) coordinates are represented in the VCF file and outputs SNPs rate per gene or a Circos heatmap track of SNP rate/bin size

### Other scripts
* `meanMedianMinMax.py`: takes input of a list of numbers and outputs the mean, median, minimum value, maximum value, and sum total
* `repeatMaskerGFFsubset`: takes input of RepeatMasker GFF and writes lines from several categories each into their own file
* `repeatMaskerGFFsummarize`: writes tables with summarized counts and lengths of features in a RepeatMasker-derived GFF3
