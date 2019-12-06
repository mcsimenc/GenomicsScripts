# Genomics scripts
Python scripts for manipulating various genomics-related file formats. \
Some definitions are included in multiple scripts so they are as standalone
as possible.

__gffAddAttribute.py__: add a key-value pair to the attributes column of a GFF3 file \
\
`gff2bed.py`: convert a GFF3 file to a BED format file \
\
`gff2circosHeatmap.py`: convert feature coordinates in a GFF3 file to Circos heatmap track format with specified bin size \
\
`gff2circosTile.py`: convert features in a GFF3 file to Circos tile track format \
\
`gff2fasta.py`: extract sequences from a FASTA file based on coordinates in a GFF3 file, using the value from a specified key in the GFF3 attributes column \
\
`gffFilter.py`: \
\
`gff2introns.py`: create a GFF3 with intron features from a GFF3 with gene and exon features or output a list of intron lengths \
\
`gffMergeOverlaps.py`: merge overlapping features in a GFF3 file \
\
`gffRemoveScafPart.py`: remove features in a GFF3 file whose coordinates \
\
`gffRenameScafs.py` \
\
