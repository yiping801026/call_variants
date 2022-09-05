dir_path=/fshare2/pingyi/cfDNA/single_call_varriants
ref=/fshare2/pingyi/ref/bwa_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bed=/fshare2/pingyi/cfDNA/NanOnCT_Panel_v1.0_Covered_hg38.bed
gatk=/fshare2/pingyi/gatk-4.2.6.1/gatk
annovar_dir=/fshare2/pingyi/bin/annovar
#mkdir $dir_path/output
cd $dir_path/try_bam
for file in `ls *duplicates.bam`
do 
    cd $dir_path/try_bam
    echo $file
    filename=$(echo $file | cut -d . -f1)
    sample_name=$(echo $filename | cut -d _ -f 1-3)
    echo $sample_name
    #1. add group to bam files：
    picard AddOrReplaceReadGroups I=$file O=$dir_path/output/"${sample_name}_group.bam" RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample_name
    #2. add.bai files for bam files：
    cd $dir_path/output/
    samtools index $dir_path/output/"${sample_name}_group.bam" 
    #3. get the middle files .gvcf
    $gatk HaplotypeCaller --emit-ref-confidence GVCF --sample-name $sample_name -R $ref -I $dir_path/output/"${sample_name}_group.bam" -L $bed -O $dir_path/output/"${sample_name}.g.vcf" && echo "** gvcf done **"
    #4. filter for gvcf
    $gatk GenotypeGVCFs \
    -R $ref \
    -V $dir_path/output/"${sample_name}.g.vcf"  \
    -L $bed \
    -O $dir_path/output/"${sample_name}.vcf"  && echo "** vcf done **"
    #5. zip 
    bgzip -f $dir_path/output/"${sample_name}.vcf"
    #6.构建tabix索引
    tabix -p vcf $dir_path/output/"${sample_name}.vcf.gz"
    #7. use SelectVariants to select SNP
    $gatk SelectVariants \
    -select-type SNP \
    -V $dir_path/output/"${sample_name}.vcf.gz" \
    -O $dir_path/output/"${sample_name}.snp.vcf.gz"

    #8. hard filter for SNPs
    $gatk VariantFiltration \
    -V $dir_path/output/"${sample_name}.snp.vcf.gz" \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O $dir_path/output/"${sample_name}.snp.filter.vcf.gz"

    #9. use SelectVariants to choose Indel
    $gatk SelectVariants \
    -select-type INDEL \
    -V $dir_path/output/"${sample_name}.vcf.gz" \
    -O $dir_path/output/"${sample_name}.indel.vcf.gz"

    #10. hard filter for Indel
    $gatk VariantFiltration \
    -V $dir_path/output/"${sample_name}.indel.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O $dir_path/output/"${sample_name}.indel.filter.vcf.gz"

    #11. bind SNP and Indel
    $gatk  MergeVcfs \
    -I $dir_path/output/"${sample_name}.snp.filter.vcf.gz" \
    -I $dir_path/output/"${sample_name}.indel.filter.vcf.gz" \
    -O $dir_path/output/"${sample_name}.filter.vcf.gz"

    #12. remove useless files 
    rm -f $dir_path/output/"${sample_name}.snp.filter.vcf.gz" \
    $dir_path/output/"${sample_name}.indel.filter.vcf.gz" \
    $dir_path/output/"${sample_name}.snp.vcf.gz" \
    $dir_path/output/"${sample_name}.indel.vcf.gz"

    #13. annotiation for vcf 
    cd annovar_dir
    perl table_annovar.pl $dir_path/output/"${sample_name}.filter.vcf" humandb/ -buildver hg38 -out "${sample_name}_final" -remove -protocol refGene,avsnp150,gnomad_genome,exac03,1000g2015aug_eas,1000g2015aug_all,clinvar_20170501,cosmic70 -operation g,f,f,f,f,f,f,f -nastring . -vcfinput -polish

    #14. move result to output_dir
    cp $annovar_dir/"${sample_name}_final.hg38_multianno.txt" $dir_path/output

    #15. filter txt to get the final variants 

done




