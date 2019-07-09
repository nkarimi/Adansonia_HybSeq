# use rename.sh to rename files
# machine : http://www.biology-it.iastate.edu/speedy2-0
# forgot to make them all R1/R2 so:
rename _ _R *

### trim files ###
# machine : http://www.biology-it.iastate.edu/speedy2-0
module load trimmomatic/0.36

for a in *_R1.fq.gz; do trimmomatic PE -threads 20 $a ${a%_R1.fq.gz}_R2.fq.gz -baseout ${a%_R1.fq.gz}.trim.fq.gz ILLUMINACLIP:Adapters.fa:2:30:15 LEADING:28 TRAILING:28 SLIDINGWINDOW:8:28 SLIDINGWINDOW:1:10 MINLEN:65 TOPHRED33; done
rename 1P R1 *
rename 2P R2 *
for a in *1U.fq.gz; do cat $a ${a%_1U.fq.gz}_2U.fq.gz > ${a%_1U.fq.gz}_S.fq.gz; done

### create bwa index and map ###
# machine : http://www.biology-it.iastate.edu/speedy2-0
module load bwa/0.7.15
bwa index GroverBaumWendelBaits_v2.fa

# split over Legion-3 (Adi*), Legion-2 (A[gmp]*), Legion-4(A[rs]*), Legion-4 (Aza*, [BCPS]*) 
# note: the Legions have only 169G total each; will have to do line-by-line for now and mv files in-between
ls *_R1.fq.gz | cut -f1 -d '.' | while read line; do bwa mem -R "@RG\\tID:$line\\tPL:ILLUMINA\\tPU:NULL\\tLB:Shotgun_LIB\\tSM:$line-NULL" -t 200 -M GroverBaumWendelBaits_v2.fa $line.trim_R1.fq.gz $line.trim_R2.fq.gz > $line.PE.sam; bwa mem -R "@RG\\tID:$line\\tPL:ILLUMINA\\tPU:NULL\\tLB:Shotgun_LIB\\tSM:$line-NULL" -t 20 -M GroverBaumWendelBaits_v2.fa $line.trim_S.fq.gz > $line.SE.sam; done

### create picard and samtools indexes required for gatk
module load picard/2.9.0
module load samtools/1.3.1.1
picard CreateSequenceDictionary R=GroverBaumWendelBaits_v2.fa O=GroverBaumWendelBaits_v2.dict
samtools faidx GroverBaumWendelBaits_v2.fa

### merge and go sam2bam
# Legion only has 169 G max HDD, do these step by step, checking output
for a in *.sam; do samtools view -bS -F 4 -q 30 $a -T GroverBaumWendelBaits_v2.fa -o ${a%.sam}.bam; samtools flagstat ${a%.sam}.bam &> ${a%.sam}.log; done

module purge
module load samtools/1.2 # avoids the error [W::bam_merge_core2] No @HD tag found
for a in *.PE.bam; do samtools merge -@ 200 ${a%.PE.bam}.bam $a ${a%.PE.bam}.SE.bam; done

for a in *.bam; do samtools sort $a ${a%.bam}.sort; done

#### purge modules and start picard/gatk ####

module load parallel/20160422
module load picard/2.9.0
module load samtools/1.2
module rm java/1.7.0_55
module load gatk/3.6

### remove duplicates, then index 
picard MarkDuplicates I=$1.sort.bam O=$1.dedup.bam REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE METRICS_FILE=$1.dedup.stats &> $1.log
samtools index $1.dedup.bam &> $1.log

ls *.bam | parallel --jobs 10 picard MarkDuplicates I={} O={.}.dedup.bam REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE METRICS_FILE={.}.dedup.stats &> {.}.log
ls *.dedup.bam | parallel samtools index {} &> {.}.log

### realign intervals
ls *.dedup.bam | parallel --jobs 10 gatk -T RealignerTargetCreator -R GroverBaumWendelBaits_v2.fa -I {} -o {.}.realign.intervals -nt 20 2> {.}.RTC.err
ls *.dedup.bam | parallel --jobs 10 gatk -T IndelRealigner -R GroverBaumWendelBaits_v2.fa -I {} -targetIntervals {.}.realign.intervals -o {.}.realign.bam 2> {.}.indel.err

### call variants and then joint genotyping

ls *.realign.bam | parallel --jobs 10 gatk -T HaplotypeCaller -R GroverBaumWendelBaits_v2.fa -I {} -nct 20 -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --genotyping_mode DISCOVERY -o {.}.raw_variants.g.vcf 2>{.}.HaplotypeCaller.err

### make sure all worked :(
mkdir done
for a in *.g.vcf; do if [[ -s $a ]]; then mv $a done; fi; done


### join the vcf for each chromosome for all accessions
# note make individual vcfs for hapcut; go to hapcut.sh
variants=""
for i in *.g.vcf; do variants+="--variant $i "; done

gatk -T GenotypeGVCFs -R GroverBaumWendelBaits_v2.fa -stand_emit_conf 10 -stand_call_conf 30 $variants --out all.Adansonia.vcf 2> joint.genotype.err


### generate hard filter for SNPs and indels #runs very fast

gatk -T SelectVariants -R GroverBaumWendelBaits_v2.fa -V all.Adansonia.vcf -selectType SNP -o all.raw_snp.vcf 2> all.SelectVariants_snp.err
 
gatk -T SelectVariants -R GroverBaumWendelBaits_v2.fa -V all.Adansonia.vcf -selectType INDEL -o all.raw_indel.vcf 2> all.SelectVariants_indel.err

# apply hard filters to snp 
 
gatk -T VariantFiltration -R GroverBaumWendelBaits_v2.fa -V all.raw_snp.vcf --filterExpression "QD<2.0||FS>60.0||MQ<40.0||SOR>4.0||MQRankSum<-12.5||ReadPosRankSum<-8.0" --filterName all_snp_filter -o all.filtered_snp.vcf 2> all.VariantFiltration_snp.err


#apply hard filters to indel 

gatk -T VariantFiltration -R GroverBaumWendelBaits_v2.fa -V all.raw_indel.vcf --filterExpression "QD<2.0||FS>200.0||SOR>10.0||ReadPosRankSum<-20.0" --filterName TANG_indel_filter -o all.filtered_indel.vcf 2> all.VariantFiltration_indel.err

 
gatk -T SelectVariants -R GroverBaumWendelBaits_v2.fa --variant all.filtered_snp.vcf --excludeFiltered -o all.PASS_snp.vcf 2> all.PASS_snp.err

gatk -T SelectVariants -R GroverBaumWendelBaits_v2.fa --variant all.filtered_indel.vcf --excludeFiltered -o all.PASS_indel.vcf 2> all.PASS_indel.err

## combine -- idk how to combine the snp/indel passed vcf intelligently yet
#gatk -T CombineVariants -V all.PASS_snp.vcf -V all.PASS_indel.vcf  -o Adansonia.snp.indel.vcf -R GroverBaumWendelBaits_v2.fa --#genotypemergeoption uniquify --filteredrecordsmergetype KEEP_IF_ALL_UNFILTERED

# filter sites
#gatk -T SelectVariants -R GroverBaumWendelBaits_v2.fa --variant Adansonia.minimal.vcf -o Adansonia.filtered.vcf --maxNOCALLfraction 0.25 
