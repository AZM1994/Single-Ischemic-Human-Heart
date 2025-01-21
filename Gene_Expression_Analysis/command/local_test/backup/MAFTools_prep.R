function prep_maftools
{
  SP_set REFERENCE_DIR="/home/yh174/reference"
  
  cat <(awk '{OFS="\t"; print $2,$3,$3,$4,$5,$1,$13,$19}' MAFTools/AD_panel.MH.loose.merged_filtered.tsv | grep -v "^Chr") <(awk '{OFS="\t"; print $2,$3,$4,$5,$6,$1,$16,$22}' MAFTools/AD_panel.indel.loose.merged_filtered.tsv | grep -v "^Chr") > MAFTools/AD_panel.both.loose.input
  table_annovar.pl MAFTools/AD_panel.both.loose.input ${REFERENCE_DIR}/human_annovar/ --buildver hg19 --outfile MAFTools/AD_panel.both.loose --otherinfo --remove --protocol refGene --operation g --nastring NA
  sed -e 's/Otherinfo/ID\tMAF\tCogdx/g' MAFTools/AD_panel.both.loose.hg19_multianno.txt | awk -F "\t" '{OFS="\t";split($6,a,";");split($9,b,";");print $1,$2,$3,$4,$5,a[1],$7,$8,b[1],$10,$11,$12,$13}' > MAFTools/AD_panel.both.loose.hg19_multianno.txt1
  mv -f MAFTools/AD_panel.both.loose.hg19_multianno.txt1 MAFTools/AD_panel.both.loose.hg19_multianno.txt
  
  awk '{OFS="\t"; print $1,$2,$2,$3,$4,$7,$6,"CHIP"}' MAFTools/blood769869-sup-document2-table2.txt | grep -v "^Chrom" > MAFTools/CHIP.hg38.input
  table_annovar.pl MAFTools/CHIP.hg38.input ${REFERENCE_DIR}/human_annovar/ --buildver hg38 --outfile MAFTools/CHIP --otherinfo --remove --protocol refGene --operation g --nastring NA
  sed -e 's/Otherinfo/ID\tMAF\tCogdx/g' MAFTools/CHIP.hg38_multianno.txt | awk -F "\t" '{OFS="\t";split($6,a,";");split($9,b,";");print $1,$2,$3,$4,$5,a[1],$7,$8,b[1],$10,$11,$12,$13}' > MAFTools/CHIP.hg38_multianno.txt1
  mv -f MAFTools/CHIP.hg38_multianno.txt1 MAFTools/CHIP.hg38_multianno.txt
}