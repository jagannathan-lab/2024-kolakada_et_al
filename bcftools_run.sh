bcftools +split-vep \
  gnomad_38_exo_stop.vcf \
  -a vep \
  -f '%CHROM	%POS	%REF	%ALT	%AF	%SYMBOL	%Gene	%Feature	%STRAND	%Codons	%Consequence	%cDNA_position	%CDS_position\n' \
  -d \
  -A tab \
  > stop_bcf.tsv
