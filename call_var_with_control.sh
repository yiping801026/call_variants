########################## README ###############################
### This script is made to call varriants from marked_duplicates.bam files 
### input files: two marked_duplicates.bam files(one tumor, one normal)
### reference files: 
### softwares： 1. samtools (get tumor.pileup & normal.pileup) 2. varscan (get .vcf files)3. annovar (annotiation)

########################## INPUT ###############################
### reference files:
ref_file=/fshare2/pingyi/ref/bwa_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
### sample files：
sample_name=NTSB0
dir_path=/fshare2/pingyi/cfDNA
tumor_bam=$dir_path/NTSB0P0001-1-MPN1_marked_duplicates.bam
normal_bam=$dir_path/NTSB0C0001-1-MPY1_marked_duplicates.bam
### annovar path:
annovar_dir=/fshare2/pingyi/bin/annovar
### bed file & MANE_txt:
bed_file="/fshare2/pingyi/cfDNA/call_varriants/NanOnCT_Panel_v1.0_Covered_hg38.bed"
MANE_txt="/fshare2/pingyi/cfDNA/call_varriants/MANE.GRCh38.v1.0.summary.txt"
### output files:
output_dir=$dir_path/call_varriants/out

########################## CODES step by step ###############################
mkdir $output_dir
### 1. samtools (get tumor.pileup & normal.pileup)
cd $output_dir
samtools mpileup -B -q 1 -f $ref_file $tumor_bam > $output_dir/"${sample_name}_tumor.pileup"
samtools mpileup -B -q 1 -f $ref_file $normal_bam > $output_dir/"${sample_name}_normal.pileup"
### 2. varscan (get .vcf files)
varscan somatic "${sample_name}_normal.pileup" "${sample_name}_tumor.pileup" "${sample_name}_varscan"  --output-vcf

### 3. annovar (annotiation)
cd $annovar_dir
perl table_annovar.pl $output_dir/"${sample_name}_varscan.snp.vcf" humandb/ -buildver hg38 -out "${sample_name}_final_snp" -remove -protocol refGene,avsnp150,gnomad_genome,exac03,1000g2015aug_eas,1000g2015aug_all,clinvar_20170501,cosmic70 -operation g,f,f,f,f,f,f,f -nastring . -vcfinput -polish
### added
perl table_annovar.pl $output_dir/"${sample_name}_varscan.indel.vcf" humandb/ -buildver hg38 -out "${sample_name}_final_indel" -remove -protocol refGene,avsnp150,gnomad_genome,exac03,1000g2015aug_eas,1000g2015aug_all,clinvar_20170501,cosmic70 -operation g,f,f,f,f,f,f,f -nastring . -vcfinput -polish


### 4. move result to output_dir
cp /fshare2/pingyi/bin/annovar/"${sample_name}_final_snp.hg38_multianno.txt" $output_dir
cp /fshare2/pingyi/bin/annovar/"${sample_name}_final_indel.hg38_multianno.txt" $output_dir
### 5.filter
### "_fianl_table_1_inbed.csv" is the result just filtered 'in bed' & 'VAF tumor bigger than normal'; 
### "_fianl_table_2_annot.csv" is the result filtered plus 'other  49:53 < 0.01' other including "gnomAD_genome_ALL","gnomAD_genome_EAS","ExAC_ALL","ExAC_EAS","X1000g2015aug_eas","X1000g2015aug_all"
### "_fianl_table_3_exononsynon.csv" is the result filtered plus exonic & nonsynonymous SNV
python /fshare2/pingyi/cfDNA/call_varriants/variants.py $output_dir/"${sample_name}_final_snp.hg38_multianno.txt" $output_dir/"${sample_name}_final_indel.hg38_multianno.txt" $bed_file $MANE_txt $output_dir/"${sample_name}_fianl_table_1_inbed.csv" $output_dir/"${sample_name}_fianl_table_2_annot.csv" $output_dir/"${sample_name}_fianl_table_3_exononsynon.csv"
