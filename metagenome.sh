
# Rawdata copy ##########################################################
workpath=myworkpath
rawdata_path=${workpath}/00.Rawdata

mkdir -p $rawdata_path
cp -r /rawdata/* ${rawdata_path}/


# 制作 filenames 
ls ${rawdata_path} | awk -F '_R' '{print $1}' | sort | uniq > ${workpath}/filenames

split -d -l 5 ${workpath}/filenames ${workpath}/filenames_

##########################################################################


# 质控 ###################################################################
workpath=myworkpath
script_path=${workpath}/Scripts
nohup ${script_path}/QualityControl.sh > fastp.log 2>&1 &

#

### QualityControl.sh ###
#!/bin/bash
tool_path=/public/home/lim/Tools
export PATH=${PATH}:${tool_path}

threads=28;

workpath=/workpath
rawdata_path=${workpath}/00.Rawdata
read_path=${workpath}/01.Cleandata

mkdir -p $read_path

for i in `cat ${workpath}/filenames`;
do

 fq_1=${Path1}/00.Rawdata/${i}_R1.fq.gz;
 fq_2=${Path1}/00.Rawdata/${i}_R2.fq.gz;
 r1=${Path2}/01.Cleandata/${i}_R1.paired.fq;
 r2=${Path2}/01.Cleandata/${i}_R2.paired.fq;

 fastp -i $fq_1 -I $fq_2 -o $r1 -O $r2 --thread $threads;

done
# 
##########################################################################



# 组装 ###################################################################
workpath=myworkpath
script_path=${workpath}/Scripts
cd $script_path
cat Assemble.sh | sed 's/myworkpath/${real_workpath}\/filename_00/g' > Assemble_00.sh
cat Assemble_00.sh | sed 's/filenames_00/filenames_01/g' > Assemble_01.sh
cat Assemble_00.sh | sed 's/filenames_00/filenames_02/g' > Assemble_02.sh
cat Assemble_00.sh | sed 's/filenames_00/filenames_03/g' > Assemble_03.sh
cat Assemble_00.sh | sed 's/filenames_00/filenames_04/g' > Assemble_04.sh
chmod +x *.sh
#
# submit job ...
#

### Assemble.sh ###
#!/bin/bash
tool_path=/public/home/lim/Tools/MEGAHIT-1.2.9-Linux-x86_64-static/bin
export PATH=${PATH}:${tool_path}

workpath=myworkpath
script_path=${workpath}/Scripts
read_path=${workpath}/01.Cleandata
asm_path=${workpath}/02.Assembly
threads=28

mkdir -p $asm_path

for i in `cat ${workpath}/filenames00`;
do

megahit -1 ${Rpath}/${i}_R1.paired.fq -2 ${Rpath}/${i}_R2.paired.fq -o ${asm_out}/${i}.asm --presets meta-large -m 0.95 --continue -t $threads

done
#
##########################################################################




# Contigs summary ###################################################################
tool_path=/public/home/lim/Tools # software : seqkit mmseq
export PATH=${PATH}:${tool_path}

contig_len_cut=500
split_num=24
threads=28

workpath=myworkpath
script_path=${workpath}/Scripts
asm_path=${workpath}/02.Assembly
gene_path=${workpath}/03.Gene

mkdir -p $asm_path $gene_path

# 重命名 contigs 文件
cd ${asm_path}
for i in `ls`; do sed "s/^>/>${i%.asm}__/" ${i}/final.contigs.fa > ${i}/final.contigs.renamed.fa ; mv ${i}/final.contigs.renamed.fa ${i}/${i%.asm}.contigs.renamed.fa; done
cd ${workpath}

# 汇总 contigs 文件
mkdir -p ${asm_path}/Summary; mv ${asm_path}/*/*.contigs.renamed.fa ${asm_path}/Summary/

# 合并 contigs
cat ${asm_path}/Summary/*.contigs.renamed.fa > ${asm_path}/contigs.renamed.fa

# filt 500 bp
seqkit seq -m ${contig_len_cut} ${asm_path}/contigs.renamed.fa > ${gene_path}/contigs.renamed.filtered.G${contig_len_cut}.fa

# 拆分 contigs 
seqkit split2 -p ${split_num} ${gene_path}/contigs.renamed.filtered.G${contig_len_cut}.fa -O ${gene_path}/contigs.renamed.filtered.G${contig_len_cut}.fa.split -f

## parallel predict ORFs
nohup find ${gene_path}/contigs.renamed.filtered.G${contig_len_cut}.fa.split/*  -name "*.fa" | xargs -n 1 -P ${split_num} ${script_path}/prodigal.sh > prodigal.log 2>&1 &

# 合并预测结果
cat ${gene_path}/contigs.renamed.filtered.G${contig_len_cut}.fa.split/*.genes.fna > ${gene_path}/contigs.renamed.filtered.${contig_len_cut}.genes.fna
cat ${gene_path}/contigs.renamed.filtered.G${contig_len_cut}.fa.split/*.proteins.faa > ${gene_path}/contigs.renamed.filtered.${contig_len_cut}.proteins.faa
cat ${gene_path}/contigs.renamed.filtered.G${contig_len_cut}.fa.split/*.genes.gbk > ${gene_path}/contigs.renamed.filtered.${contig_len_cut}.genes.gbk 


### prodigal.sh ###
#!/bin/bash
export PATH=${PATH}:/public/home/lim/Tools/Prodigal-v2.6.3

input_file="$1"

prodigal -i $input_file -a ${input_file%.*}.proteins.faa -d ${input_file%.*}.genes.fna -o ${input_file%.*}.genes.gbk -p meta -q 
#

###################################################################





# 基因聚类 ###################################################################
workpath=myworkpath
script_path=${workpath}/Scripts
cd $script_path
cat cluster.sh | sed "s|myworkpath|$workpath|g" > cluster_.sh
chmod +x cluster_.sh
nohup ./cluster_.sh > mmseq.log 2>&1 &
#

### cluster.sh ###
#!/bin/bash
source /public/home/lim/Tools/anaconda3/etc/profile.d/conda.sh
conda activate mmseq

contig_len_cut=500
threads=28

workpath=myworkpath
gene_path=${workpath}/03.Gene

cd $gene_path

mkdir -p DB DB_90

mmseqs createdb contigs.renamed.filtered.G${contig_len_cut}.proteins.faa DB/DB

mmseqs cluster --threads ${threads} --min-seq-id 0.90 -c 0.90 --cov-mode 1 --cluster-mode 2 DB/DB DB_90/DB_90 tmp_90

mmseqs createtsv --threads ${threads} DB/DB DB/DB DB_90/DB_90 DB_90.tsv

##########################################################################




# extract representative sequences (pipeline_05.sh) #############################################
contig_len_cut=500
threads=28
workpath=myworkpath
gene_path=${workpath}/03.Gene

cd $gene_path

cut -f 1 DB_90.tsv | sort | uniq > DB_90.representative.tsv

seqkit replace -p " .*" -r "" contigs.renamed.filtered.G${contig_len_cut}.genes.fna > contigs.renamed.filtered.G${contig_len_cut}.genes.cleaned.fna
seqkit replace -p " .*" -r "" contigs.renamed.filtered.G${contig_len_cut}.proteins.faa > contigs.renamed.filtered.G${contig_len_cut}.proteins.cleaned.faa

awk 'NR==FNR {id[$1]; next} /^>/ {f=0} substr($1,2) in id {f=1} f' DB_90.representative.tsv contigs.renamed.filtered.G${contig_len_cut}.genes.cleaned.fna > contigs.renamed.filtered.G${contig_len_cut}.genes.cleaned.representative.fna
awk 'NR==FNR {id[$1]; next} /^>/ {f=0} substr($1,2) in id {f=1} f' DB_90.representative.tsv contigs.renamed.filtered.G${contig_len_cut}.proteins.cleaned.faa > contigs.renamed.filtered.G${contig_len_cut}.proteins.cleaned.representative.faa

# optional #
seqkit seq -j ${threads} -m 150 contigs.renamed.filtered.G${contig_len_cut}.genes.cleaned.representative.fna > contigs.renamed.filtered.G${contig_len_cut}.genes.cleaned.representative.cut150.fna
seqkit seq -j ${threads} -m 50 contigs.renamed.filtered.G${contig_len_cut}.proteins.cleaned.representative.faa > contigs.renamed.filtered.G${contig_len_cut}.proteins.cleaned.representative.cut50.faa

################################################################################





# Gene KO annotation ######################################################
workpath=myworkpath
gene_path=${workpath}/03.Gene
script_path=${workpath}/Scripts

cd $gene_path
mkdir -p ${gene_path}/KO.scan
cat ${gene_path}/contigs.renamed.filtered.G500.proteins.cleaned.representative.cut50.faa | seqkit split2 -s 100000 -O ${gene_path}/KO.scan -f

cd $script_path
vim koscan.sh
cat koscan.sh | sed "s|myworkpath|$workpath|g" > koscan_.sh
cat koscan_.sh | sed 's/stdin.part_/stdin.part_0/g' > koscan_00.sh
cat koscan_.sh | sed 's/stdin.part_/stdin.part_1/g' > koscan_01.sh
cat koscan_.sh | sed 's/stdin.part_/stdin.part_2/g' > koscan_02.sh
cat koscan_.sh | sed 's/stdin.part_/stdin.part_3/g' > koscan_03.sh
cat koscan_.sh | sed 's/stdin.part_/stdin.part_4/g' > koscan_04.sh
cat koscan_.sh | sed 's/stdin.part_/stdin.part_5/g' > koscan_05.sh
cat koscan_.sh | sed 's/stdin.part_/stdin.part_6/g' > koscan_06.sh
chmod +x koscan_*.sh
#
# submit job ...
#

# koscan.sh
#!/bin/bash
source /public/home/lim/Tools/anaconda3/etc/profile.d/conda.sh 
conda activate kofam_scan

export PATH=$PATH:/public/home/lim/Tools/kofam_scan-1.3.0

workpath=myworkpath
gene_path=${workpath}/03.Gene

for p in ${gene_path}/KO.scan/stdin.part_*;
do

n="${p%.*}" ;

exec_annotation -p /public/home/lim/Database/kegg/profiles -k /public/home/lim/Database/kegg/ko_list -f detail-tsv --tmp-dir ${n}.tmp -o ${n}.queryKOs $p

done
#
################################################################################



# koscan result clean ##########################################################
workpath=myworkpath
gene_path=${workpath}/03.Gene
script_path=${workpath}/Scripts
cd $script_path
vim clean.sh; chmod +x clean.sh 
find ${gene_path}/KO.scan -name "*.queryKOs" | xargs -n 1 -P 28 ./clean.sh
cat ${gene_path}/KO.scan/*_processed.txt > ${gene_path}/queryKOs.txt
awk '{ if($6 < 0.00005) {print $0} }' ${gene_path}/queryKOs.txt > ${gene_path}/queryKOs_filt.txt 
#

### clean.sh ###
# 处理 kofam_scan 输出结果
input_file="$1"
output_file="${input_file}_processed.txt"

grep "^*" "$input_file" > "$output_file"
#
################################################################################






# Salmon gene quatify ########################################################################
contig_len_cut=500
scp lim@10.75.73.96:${workpath}/03.Gene/contigs.renamed.filtered.G${contig_len_cut}.genes.cleaned.representative.cut150.fna ${workpath}/03.Gene/
scp -r lim@10.75.73.96:${workpath}/01.Cleandata ${workpath}/01.Cleandata

workpath=myworkpath
script_path=${workpath}/Scripts
cd $script_path
vim salmon_index.sh
cat salmon_index.sh | sed "s|myworkpath|$workpath|g" > salmon_index_.sh

vim salmon.sh
cat salmon.sh | sed "s|myworkpath|$workpath|g" > salmon_.sh
cat salmon_.sh | sed 's/filenames/filenames_00/g' > salmon_00.sh
cat salmon_.sh | sed 's/filenames/filenames_01/g' > salmon_01.sh
cat salmon_.sh | sed 's/filenames/filenames_02/g' > salmon_02.sh
cat salmon_.sh | sed 's/filenames/filenames_03/g' > salmon_03.sh
cat salmon_.sh | sed 's/filenames/filenames_04/g' > salmon_04.sh
chmod +x *.sh
#
# submit job ...
#

### salmon_index.sh ###
#!/bin/bash

source /public/home/lim/micromamba/etc/profile.d/micromamba.sh
micromamba activate salmon

workpath=myworkpath
gene_path=${workpath}/03.Gene

salmon index -t ${gene_path}/contigs.renamed.filtered.G500.genes.cleaned.representative.cut150.fna -i ${gene_path}/salmon_index -p 64
#

### salmon.sh ###
#!/bin/bash

source /public/home/lim/micromamba/etc/profile.d/micromamba.sh
micromamba activate salmon

workpath=myworkpath
read_path=${workpath}/01.Cleandata
gene_path=${workpath}/03.Gene

for id in `cat ${workpath}/filenames`
do

salmon quant -i ${gene_path}/salmon_index -l A -1 ${read_path}/${id}_R1.paired.fq -2 ${read_path}/${id}_R2.paired.fq -p 64 --meta -o ${gene_path}/${id}.quant

done
#
###############################################################################################




############### merge salmon results ######################################################
workpath=myworkpath
gene_path=${workpath}/03.Gene
cd $gene_path

mkdir -p temp_output

# 遍历当前目录下所有子文件夹中的 result.txt
for dir in *.quant/ ; do
    file="${dir}/quant.sf"
    if [[ -f "$file" ]]; then
        # 提取倒数第二列，保存到临时文件中
        awk '{print $(NF-1)}' "$file" > "temp_output/${dir%/}.txt"
    fi
done

cd temp_output

for file in *.txt; do
    filename=${file%%.*}
    sed -i "1s/.*/$filename/" "$file"
done

cd ..

cut -f1 "$(ls */quant.sf | head -n1)" > temp_output/00.names.txt

paste temp_output/*.txt > merged.txt

################################################################################################




## merge table and annotation ##################################################################

workpath=myworkpath
gene_path=${workpath}/03.Gene
cd $gene_path
# scp -r lim@10.75.73.96:$gene_path/queryKOs_filt.txt lim@10.75.73.77:$gene_path/
micromamba activate r-base-4.3
R

#### R ###
library(tidyverse)

table = read.csv("merged.txt", sep = '\t', header = TRUE)

queryKOs = read.csv("queryKOs_filt.txt", sep = '\t', header = FALSE)

colnames(queryKOs)[2] = "Name"

df_new = table %>% left_join(queryKOs, by = "Name")

cols_to_sum <- names(df_new)[ !startsWith(names(df_new), "V") & names(df_new) != "Name" ]

df_sum <- df_new %>%
  group_by(V3) %>%
  summarise(desc = first(V7), across(all_of(cols_to_sum), sum), .groups = "drop")

head(df_sum)
ncol(df_sum)
nrow(df_sum)

colnames(df_sum)[1] = "KO"
write.csv(df_sum, "kegg_gene.csv", row.names = FALSE, quote = FALSE)

#### R ###
################################################################################################






## merge Other table and annotation ##################################################################

workpath=myworkpath
gene_path=${workpath}/03.Gene
cd $gene_path
# scp -r lim@10.75.73.96:$gene_path/queryKOs_filt.txt lim@10.75.73.77:$gene_path/
micromamba activate r-base-4.3
R

#### R ###
library(tidyverse)

table = read.csv("/public/home/lim/Dataset/zmdai20250304/TempFile/gene.txt", sep = '\t', header = TRUE)

queryKOs = read.csv("/public/home/lim/Dataset/zmdai20250304/07.MGEs_MobileOGdb/result_res", sep = '\t', header = FALSE)

colnames(queryKOs)[1] = "Name"

df_new = table %>% left_join(queryKOs, by = "Name")

cols_to_sum <- names(df_new)[ !startsWith(names(df_new), "V") & names(df_new) != "Name" ]

df_sum <- df_new %>%
  group_by(V2) %>%
  summarise(Evalue = first(V11), across(all_of(cols_to_sum), sum), .groups = "drop")

head(df_sum)
ncol(df_sum)
nrow(df_sum)

colnames(df_sum)[1] = "Gene"
write.csv(df_sum, "MGEs_gene.csv", row.names = FALSE, quote = FALSE)

#### R ###
################################################################################################







# 单样本 Mapping ###############################################################
workpath=myworkpath
script_path=${workpath}/Scripts
cd $script_path
vim bowtie.mapping.sh
cat bowtie.mapping.sh | sed "s|myworkpath|$workpath|g" > bowtie.mapping_.sh 
cat bowtie.mapping_.sh | sed 's/filenames/filenames_00/g' > bowtie.mapping_01.sh
cat bowtie.mapping_.sh | sed 's/filenames/filenames_01/g' > bowtie.mapping_02.sh
cat bowtie.mapping_.sh | sed 's/filenames/filenames_02/g' > bowtie.mapping_03.sh
cat bowtie.mapping_.sh | sed 's/filenames/filenames_03/g' > bowtie.mapping_04.sh
cat bowtie.mapping_.sh | sed 's/filenames/filenames_04/g' > bowtie.mapping_05.sh
chmod +x *.sh
#
# submit job ...
#

### bowtie.mapping.sh ###
#!/bin/bash
export PATH=/public/home/lim/Tools/bowtie2-2.5.1-linux-x86_64:${PATH}
export PATH=/public/home/lim/Tools/samtools-1.15.1:${PATH}

workpath=myworkpath
read_path=${workpath}/01.Cleandata
asm_path=${workpath}/02.Assembly
bin_path=${workpath}/04.Binning

for n in `cat ${workpath}/filenames`
do 
	if [ -f "${read_path}/${n}_R1.paired.fq" ]
	then
		mkdir -p ${bin_path}/Align/${n}/index

		bowtie2-build ${asm_path}/${n}.contigs.renamed.filtered.fa  ${bin_path}/Align/${n}/index/${n}.index

  		bowtie2 -x ${bin_path}/Align/${n}/index/${n}.index -1 ${read_path}/${n}_R1.paired.fq -2 ${read_path}/${n}_R2.paired.fq -p 64 -S ${bin_path}/Align/${n}/${n}.sam 

		samtools view -h -b -F 4 -S ${bin_path}/Align/${n}/${n}.sam -o ${bin_path}/Align/${n}/${n}.bam

		rm ${bin_path}/Align/${n}/${n}.sam
		
		samtools sort -@ 64 ${bin_path}/Align/${n}/${n}.bam -o ${bin_path}/Align/${n}/${n}.sorted.bam
	fi
done
#
###############################################################################








# 单样本 semibin binning ##############################################################
workpath=myworkpath
script_path=${workpath}/Scripts
cd $script_path
vim semibin.pipeline.sh
cat semibin.pipeline.sh | sed "s|myworkpath|$workpath|g" > semibin.pipeline_.sh 
cat semibin.pipeline_.sh | sed 's/filenames/filenames_00/g' > semibin.pipeline_00.sh
cat semibin.pipeline_.sh | sed 's/filenames/filenames_01/g' > semibin.pipeline_01.sh
cat semibin.pipeline_.sh | sed 's/filenames/filenames_02/g' > semibin.pipeline_02.sh
cat semibin.pipeline_.sh | sed 's/filenames/filenames_03/g' > semibin.pipeline_03.sh
cat semibin.pipeline_.sh | sed 's/filenames/filenames_04/g' > semibin.pipeline_04.sh
chmod +x *.sh
#
# submit job ...
#

# semibin.pipeline.sh ##############################################################
#!/bin/bash
source /public/home/lim/Tools/anaconda3/etc/profile.d/conda.sh

conda activate semibin-1.5.1

workpath=myworkpath
asm_path=${workpath}/02.Assembly
bin_path=${workpath}/04.Binning

for sample in `cat ${workpath}/filenames`
do

mkdir -p ${bin_path}/advanced_single_sample_output/${sample}

# Generate data.csv/data_split.csv
SemiBin generate_sequence_features_single \
-i ${asm_path}/${sample}.contigs.renamed.filtered.fa \
-b ${bin_path}/Align/${sample}/${sample}.sorted.bam \
-o ${bin_path}/advanced_single_sample_output/${sample} \
-p 24

# Generate cannot-link
SemiBin generate_cannot_links \
-r /public/home/lim/Database/semibin-database \
-i ${asm_path}/${sample}.contigs.renamed.filtered.fa  \
-o ${bin_path}/advanced_single_sample_output/${sample} \
-p 24

# Bin 2
SemiBin bin \
-i ${asm_path}/${sample}.contigs.renamed.filtered.fa  \
--data ${bin_path}/advanced_single_sample_output/${sample}/data.csv \
-o ${bin_path}/advanced_single_sample_output_pretrain.bins/${sample} \
--environment soil \
-p 24

done
#
###########################################################################





# semibin bins merge ######################################################
workpath=myworkpath
bin_path=${workpath}/04.Binning

cd ${bin_path}/SemiBin_work

mkdir -p ./merged

find ./*/output_recluster_bins/* -type f | while read file; do
    rel_path="${file#./}"
    top_dir=$(echo "$rel_path" | cut -d'/' -f1)
    base_name=$(basename "$file")
    new_name="${top_dir}_${base_name}.semibin"
    cp "$file" "./merged/$new_name"
done

###########################################################################



# merge bins from multiple methods ########################################
workpath=myworkpath
bin_path=${workpath}/04.Binning

mkdir -p ${bin_path}/merged

cp  ${bin_path}/*_work/merged/*.fa  ${bin_path}/merged/

###########################################################################



# 去冗余 ###################################################################
workpath=myworkpath
script_path=${workpath}/Scripts; cd $script_path
vim drep.sh
cat drep.sh | sed "s|myworkpath|$workpath|g" > drep_.sh; chmod +x drep_.sh
nohup ./drep_.sh > drep.log 2>&1 &
#

### drep ###
#!/bin/bash
source /public/home/lim/micromamba/etc/profile.d/micromamba.sh
micromamba activate dRep-3.4

export CHECKM_DATA_PATH=/public/home/lim/Database/checkm_data_2015_01_16
export PATH=$PATH:/usr/local/bin

threads=28

workpath=myworkpath
bin_path=${workpath}/04.Binning

dRep dereplicate ${bin_path}/derep_out -p ${threads} -comp 50 -con 10 -g ${bin_path}/merged/*.fa

#
###########################################################################




# prodigal predict ORFs from each MAG ######################################
workpath=myworkpath
bin_path=${workpath}/04.Binning
script_path=${workpath}/Scripts; cd $script_path

vim prodigal.sh; chmod +x prodigal.sh 

Parallel=28
nohup find ${bin_path}/derep_out/dereplicated_genomes/* -name "*.fa" | xargs -n 1 -P ${Parallel} ./prodigal.sh > prodigal.log 2>&1 &

### prodigal.sh ###
#!/bin/bash
export PATH=${PATH}:/public/home/lim/Tools/Prodigal-v2.6.3

input_file="$1"

echo "Processing $input_file"

prodigal -i $input_file -o ${input_file%.*}.gene.gbk -a ${input_file%.*}.prot.faa -d ${input_file%.*}.gene.fna -p single -q 
# 
###########################################################################



# ko annotate for each MAG ######################################
workpath=myworkpath
bin_path=${workpath}/04.Binning
script_path=${workpath}/Scripts; cd $script_path

Parallel=40

find ${bin_path}/derep_out/dereplicated_genomes/*.prot.faa -type f > ${script_path}/all_files.txt
split -n l/${Parallel} ${script_path}/all_files.txt ${script_path}/files_batch_ -d

vim koscan.sh
cat koscan.sh | sed "s|myworkpath|$workpath|g" > koscan_.sh
cat koscan_.sh | sed 's/files_batch/files_batch_00/g' > koscan_00.sh
cat koscan_.sh | sed 's/files_batch/files_batch_01/g' > koscan_01.sh
cat koscan_.sh | sed 's/files_batch/files_batch_02/g' > koscan_02.sh
cat koscan_.sh | sed 's/files_batch/files_batch_03/g' > koscan_03.sh
cat koscan_.sh | sed 's/files_batch/files_batch_04/g' > koscan_04.sh
cat koscan_.sh | sed 's/files_batch/files_batch_05/g' > koscan_05.sh
cat koscan_.sh | sed 's/files_batch/files_batch_06/g' > koscan_06.sh
cat koscan_.sh | sed 's/files_batch/files_batch_07/g' > koscan_07.sh
cat koscan_.sh | sed 's/files_batch/files_batch_08/g' > koscan_08.sh
cat koscan_.sh | sed 's/files_batch/files_batch_09/g' > koscan_09.sh
cat koscan_.sh | sed 's/files_batch/files_batch_10/g' > koscan_10.sh
cat koscan_.sh | sed 's/files_batch/files_batch_11/g' > koscan_11.sh
cat koscan_.sh | sed 's/files_batch/files_batch_12/g' > koscan_12.sh
cat koscan_.sh | sed 's/files_batch/files_batch_13/g' > koscan_13.sh
cat koscan_.sh | sed 's/files_batch/files_batch_14/g' > koscan_14.sh
cat koscan_.sh | sed 's/files_batch/files_batch_15/g' > koscan_15.sh
cat koscan_.sh | sed 's/files_batch/files_batch_16/g' > koscan_16.sh
cat koscan_.sh | sed 's/files_batch/files_batch_17/g' > koscan_17.sh
cat koscan_.sh | sed 's/files_batch/files_batch_18/g' > koscan_18.sh
cat koscan_.sh | sed 's/files_batch/files_batch_19/g' > koscan_19.sh
chmod +x koscan*.sh
#
# submit job ...
#

### koscan.sh ###
#!/bin/bash
source /public/home/lim/Tools/anaconda3/etc/profile.d/conda.sh # source /public/home/lim/micromamba/etc/profile.d/micromamba.sh
conda activate kofam_scan # micromamba activate kofam_env
export PATH=$PATH:/public/home/lim/Tools/kofam_scan

workpath=myworkpath
script_path=${workpath}/Scripts

for g in `cat ${script_path}/files_batch`
do

exec_annotation -f detail-tsv -p /public/home/lim/Database/kegg/profiles -k /public/home/lim/Database/kegg/ko_list --tmp-dir ${g%.*}.tmp -o ${g%.*}.queryKOs ${g}

done
#
###########################################################################



# koscan result clean ##########################################################
workpath=myworkpath
bin_path=${workpath}/04.Binning
script_path=${workpath}/Scripts; cd $script_path

vim clean.sh; chmod +x clean.sh 

find ${bin_path}/KO.scan -name "*.queryKOs" | xargs -n 1 -P 28 ./clean.sh


#

### clean.sh ###
# 处理 kofam_scan 输出结果
input_file="$1"
output_file="${input_file}_filt.txt"

awk '$1 ~ /^\*/ && $6 <= 0.00005' "$input_file" > "$output_file"
################################################################################



# KO Matrix ########### (optional) ############################################


### R ###

library(data.table)

# 文件名列表（替换成你的文件路径）
setwd("")

files <- list.files(pattern = "queryKOs_filt\\.txt$")

all_kos <- c()

# 遍历所有文件
for (f in files) {
  # 读取文件（假设为制表符分隔）
  df <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # 提取第三列并添加到列表
  if (ncol(df) >= 3) {
    all_kos <- c(all_kos, df[[3]])
  }
}

# 去重
all_kos <- unique(sort(all_kos))

# 创建空矩阵，行数=文件数，列数=全集KO数
mat <- matrix(0, nrow=length(files), ncol=length(all_kos))
rownames(mat) <- files
colnames(mat) <- all_kos

# 填充矩阵
for(i in seq_along(files)) {
  ko <- fread(files[i], header=FALSE, sep="\t", data.table=FALSE)[,3]
  ko <- ko[ko != ""]  # 去空
  mat[i, all_kos %in% ko] <- 1
}

# 查看结果
head(mat)
row.names(mat) = sub("_[^_]*$", "", row.names(mat))

# 可选：写出到文件
write.table(mat, file="KO_matrix.tsv", sep="\t", quote=FALSE, col.names=NA)
#
###########################################################################




# CheckM2 ####################################################
workpath=myworkpath
bin_path=${workpath}/04.Binning
script_path=${workpath}/Scripts; cd $script_path
vim checkm2.sh
cat checkm2.sh | sed "s|myworkpath|$workpath|g" > checkm2_.sh

### checkm2.sh ###
#!/bin/bash
source /public/home/lim/micromamba/etc/profile.d/micromamba.sh
micromamba activate checkm2-v1.0.2

workpath=myworkpath
bin_path=${workpath}/04.Binning
threads=28

mkdir -p ${bin_path}/checkm2

checkm2 predict --threads $threads --input ${bin_path}/derep_out/dereplicated_genomes/ \
--output-directory ${bin_path}/checkm2 \
--database_path /public/home/lim/Database/checkm2/CheckM2_database/uniref100.KO.1.dmnd -x fa
#
###########################################################################





# MAG Taxonomy ####################################################
workpath=myworkpath
bin_path=${workpath}/04.Binning
script_path=${workpath}/Scripts; cd $script_path
vim gtdbtk.sh
cat gtdbtk.sh | sed "s|myworkpath|$workpath|g" > gtdbtk_.sh
chmod +x gtdbtk_.sh

### gtdbtk.sh ###
#!/bin/bash
source /public/home/lim/micromamba/etc/profile.d/micromamba.sh
micromamba activate gtdbtk-2.4.0

workpath=myworkpath
bin_path=${workpath}/04.Binning
threads=64

mkdir -p ${bin_path}/Taxonomy

gtdbtk classify_wf --genome_dir ${bin_path}/derep_out/dereplicated_genomes --out_dir ${bin_path}/Taxonomy \
--cpus $threads --extension fa \
--mash_db /public/home/lim/micromamba/envs/gtdbtk-2.3.2/share/gtdbtk-2.3.2/db/release214/mash/mash.msh
#
###########################################################################




# MAG Abundance #####################################################
workpath=myworkpath
bin_path=${workpath}/04.Binning
script_path=${workpath}/Scripts; cd $script_path

read_path=${workpath}/01.Cleandata
ls $read_path | awk -v path="\${read_path}" '{print path "/" $0}' | tr '\n' ' '
vim coverm.sh
cat coverm.sh | sed "s|myworkpath|$workpath|g" > coverm_.sh

### coverm.sh ###
#!/bin/bash
source /public/home/lim/micromamba/etc/profile.d/micromamba.sh
micromamba activate coverm

workpath=myworkpath
read_path=${workpath}/01.Cleandata
bin_path=${workpath}/04.Binning
threads=28

coverm genome -d ${bin_path}/derep_out/dereplicated_genomes -x fa -t $threads --methods relative_abundance \
-c ${read_path}/*_R1.paired.fq ${read_path}/*_R2.paired.fq \
-o ${bin_path}/bin_relative_abundance.tsv

#% MAG covered
coverm genome \
  --proper-pairs-only \
  --genome-fasta-extension fa \
  --genome-fasta-directory ${bin_path}/derep_out/dereplicated_genomes  \
  -c ${read_path}/*_R1.paired.fq ${read_path}/*_R2.paired.fq \
  --methods rpkm tpm relative_abundance covered_fraction trimmed_mean reads_per_base count \
  --threads 64 \
  --min-read-percent-identity-pair 0.95 \
  --min-covered-fraction 0.75 \
  --output-file ${bin_path}/bowtie_mapping95id_MAGs_MAGcovered_output
#
###########################################################################






