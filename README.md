# Genomics scripts
Python scripts for manipulating various genomics-related file formats. \
Some definitions are included in multiple scripts so they are as standalone
as possible.

### BED scripts
* `bedIntersect2percentOverlap.py`: output proportion of features participating in an overlap as found by bedtools bedIntersect -wa -wb
* `gff2bed.py`: convert a GFF3 file to a BED format file

### BLAST scripts
* `blast2gff`: convert blastn, blastp, etc. tabular output to GFF3 format
* `blastBestHit.py`: outputs the highest-scoring hit from blastn, blastp, etc. tabular output (-outfmt 6 or 7)
* `blastFilter.py`: outputs lines where the percent identity satisfies constraints from blastn, blastp, etc. tabular output (-outfmt 6 or 7)

### Circos scripts
* `coverage2circosLine.py`: calculate average depth of coverage form the output of bedtools genomecov -ibam <bam> -d and output a Circos line track
* `gff2circosHeatmap.py`: convert feature coordinates in a GFF3 file to Circos heatmap track format with specified bin size
* `gff2circosTile.py`: convert features in a GFF3 file to Circos tile track format

### FASTA scripts
* `extractFastaSeqs.py`: extract a subset of the sequences in a FASTA file
* `gff2fasta.py`: extract sequences from a FASTA file based on coordinates in a GFF3 file, using the value from a specified key in the GFF3 attributes column

### GFF scripts
* `blast2gff`: convert blastn, blastp, etc. tabular output to GFF3 format
* `gff2bed.py`: convert a GFF3 file to a BED format file
* `gff2circosHeatmap.py`: convert feature coordinates in a GFF3 file to Circos heatmap track format with specified bin size
* `gff2circosTile.py`: convert features in a GFF3 file to Circos tile track format
* `gff2fasta.py`: extract sequences from a FASTA file based on coordinates in a GFF3 file using the value from a specified key in the GFF3 attributes column as the sequence name. Depends on BEDTools and BioPython
* `gff2introns.py`: create a GFF3 with intron features from a GFF3 with gene and exon features or output a list of intron lengths
* `gffAddAttribute.py`: add a key-value pair to the attributes column of a GFF3 file
* `gffFilter.py`: remove or retain GFF3 features on specified scaffolds or with specified values for the ID attribute
* `gffMergeOverlaps.py`: merge overlapping features in a GFF3 file
* `gffRemoveScafPart.py`: remove features in a GFF3 file whose coordinates
* `gffRenameScafs.py`: rename scaffolds in a GFF3 file per a two-column map
