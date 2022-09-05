import sys 
import pandas as pd
from pandas import DataFrame

def filter_column(final_txt,col_orig,col_filter):
    for n in range(1,len(final_txt)):
        if final_txt[col_orig][n]=='.':
            final_txt[col_filter][n] = 'TRUE'
        else:
            if float(final_txt[col_orig][n]) < 0.01:
                final_txt[col_filter][n] = 'TRUE'
            else:
                final_txt[col_filter][n] = 'FALSE'
    return final_txt


def filter_base (bed_file,final_txt):
    ### 49 for bed 
    final_txt[49]=final_txt[0]
    ### 50:53 for threshold < 0.01
    final_txt[50]=final_txt[0]
    final_txt[51]=final_txt[0]
    final_txt[52]=final_txt[0]
    final_txt[53]=final_txt[0]
    final_txt[54]=final_txt[0]
    
    final_txt = filter_column(final_txt,11,49)
    final_txt = filter_column(final_txt,15,50)
    final_txt = filter_column(final_txt,19,51)
    final_txt = filter_column(final_txt,22,52)
    final_txt = filter_column(final_txt,27,53)
    final_txt = filter_column(final_txt,28,54)
    ### in bed or not ###
    final_txt[55]=final_txt[0]
    for i in range(1,len(final_txt)):
        for m in range(len(bed_file)):
            if final_txt[0][i] == bed_file[0][m] and int(final_txt[1][i]) in range(bed_file[1][m],bed_file[2][m]+1) and final_txt[6][i] == bed_file[3][m]:
                final_txt[55][i] = '1'
    
    ### filter in bed(55) & exonic & nonsynonymous SNV & 49:53 < 0.01
    ### only filter in bed 
    final_txt = final_txt.loc[final_txt[55].isin(['1'])]

    return final_txt

def get_transcript (final_filter,MANE_txt_table):
    rownames = final_filter.index
    ###
    final_filter[56] = final_filter[0]
    final_filter[57] = final_filter[0]
    final_filter[58] = final_filter[0]
    final_filter[59] = final_filter[0]
    final_filter[60] = final_filter[0]
    for p in rownames:
        RefSeq = MANE_txt_table.loc[MANE_txt_table['symbol'].isin([final_filter[6][p]]),'RefSeq_nuc']
        RefSeq = str(RefSeq)
        final_filter[56][p] = RefSeq.split('\n')[0].split(' ')[-1].split('.')[0]
        ### choose transcript
        list_seq = final_filter[9][p].split(',')
        
        for l in list_seq: 
            if l == '.':
                final_filter[59][p] = '.'
            else:
                if final_filter[56][p] == l.split(':')[1]:
                    final_filter[59][p] = l
        ### 47 < 48 VAF
        if float(final_filter[47][p].split(':')[5].strip("%"))<float(final_filter[48][p].split(':')[5].strip("%")):
            final_filter[57][p] = 'TRUE'
        ### 58 is the VAF of tumor
        final_filter[58][p] = final_filter[48][p].split(':')[5]

        ### 60 is germline <= 30%
        if float(final_filter[47][p].split(':')[5].strip("%")) < 30:
            final_filter[60][p] = 'TRUE'
        else:
            final_filter[60][p] = 'FALSE'

    final_filter = final_filter.loc[final_filter[57].isin(['TRUE'])]
    return final_filter



#### ------ main ------ ####
#### ------ read tables ------ ####
### read sample
final_snp_txt = pd.read_table(sys.argv[1],sep = '\t',header=None)
final_indel_txt = pd.read_table(sys.argv[2],sep = '\t',header=None)

### merge snp & indel together
final_txt = pd.concat([final_snp_txt,final_indel_txt[1:]])

### change the row names for table
final_txt = pd.DataFrame(final_txt)
final_txt.index=range(0,len(final_txt))

### read other
bed_file = pd.read_table(sys.argv[3],sep = '\t',header=None)
MANE_txt = pd.read_table(sys.argv[4],sep = '\t')

####### filter just 1. in bed && 2. VAF tumor > normal ,
### output: '~/Amazing_Yi/work/TH/projects/callvarriants/variants/fianl_table004_1_bedVAF.csv'
final_base_filter = filter_base(bed_file,final_txt)

### choose the VAF & transcript(AAChange.refGene in ) & choose tumor bigger than normal '48>47'
final_trans =  get_transcript(final_base_filter,MANE_txt)

####### filter + with other  49:54 < 0.01, 
fianl_table004_2_annot = final_trans.loc[final_trans[49].isin(['TRUE']) & final_trans[50].isin(['TRUE']) & final_trans[51].isin(['TRUE']) & final_trans[52].isin(['TRUE']) & final_trans[53].isin(['TRUE'])& final_trans[54].isin(['TRUE'])]
####### filter + exonic & nonsynonymous SNV
fianl_table004_3_exononsynon = fianl_table004_2_annot.loc[fianl_table004_2_annot[5].isin(['exonic']) & fianl_table004_2_annot[8].isin(['nonsynonymous SNV'])]

############## print results
col_list = [0,1,2,3,4,5,6,7,8,
44,45,46,47,48,49,50,51,52,53,54,57,
58,59,
60]

colname_list = ["Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene",
"Qulit","Otherinfo11","Otherinfo12","Otherinfo13","Otherinfo14","gnomAD_genome_ALL","gnomAD_genome_EAS","ExAC_ALL","ExAC_EAS","X1000g2015aug_eas","X1000g2015aug_all","VAF_tumor_bigger",
"VAF","transcript",
"germ_samller30%" ]

######## only in bed + tumor bigger
final_table = final_trans.iloc[:,col_list]
final_table.columns = colname_list
final_table.to_csv(sys.argv[5])

####### filter + with other  49:53 < 0.01, 
fianl_table004_2_annot = fianl_table004_2_annot.iloc[:,col_list]
fianl_table004_2_annot.columns = colname_list
fianl_table004_2_annot.to_csv(sys.argv[6])

####### filter + exonic & nonsynonymous SNV
fianl_table004_3_exononsynon = fianl_table004_3_exononsynon.iloc[:,col_list]
fianl_table004_3_exononsynon.columns = colname_list
fianl_table004_3_exononsynon.to_csv(sys.argv[7])


