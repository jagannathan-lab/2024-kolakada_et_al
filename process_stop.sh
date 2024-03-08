#sed  '/>/ s/\..*//' /Users/rf/Downloads/gencode.v19.pc_transcripts.fa > gencode.v19.pc_transcripts_renamed.fa 
    #fix names for seq lookup
#samtools faidx gencode.v19.pc_transcripts_renamed.fa 
    #index

awk '$16 ~ /stop/' chr1_stop.tsv | awk 'length($3) ==1 && length($4) ==1' | awk '$3 != "-" && $4 != "-"' > chr1_stop_filtered.tsv 
    #only keeping single nucleotide variants with stop_gained
wc -l chr1_stop_filtered.tsv

#constructing bed format, and retrieve original codon, minus1, and plus1
awk -v OFS='\t' '{$32 = $17 - ($18 + 2) % 3; $33 = $32 + 2}1' chr1_stop_filtered.tsv | \
  awk '$32 >=4' | \
  awk -v OFS='\t' '{print $13, $32 - 1, $33, $0}' | \
  bedtools getfasta -fi gencode.v19.pc_transcripts_renamed.fa -bed - -bedOut | \
  awk -v OFS='\t' '{print $1, $2 - 3, $3 - 3, $0}' | \
  bedtools getfasta -fi gencode.v19.pc_transcripts_renamed.fa -bed - -bedOut | \
  awk -v OFS='\t' '{print $1, $2 + 6, $3 + 6, $0}' | \
  bedtools getfasta -fi gencode.v19.pc_transcripts_renamed.fa -bed - -bedOut > chr1_stop_codons.tsv
wc -l chr1_stop_codons.tsv
