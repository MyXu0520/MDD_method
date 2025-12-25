zcat gencode.v44.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sort -k1,1 -k2,2n| bedtools merge > hg38.exon.sort.bedtools.merge.bed
# 1.
zcat gencode.v44.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sort -k1,1V -k2,2n >hg38.gene.sort.bed
# 2.
bedtools subtract -a hg38.gene.sort.bed -b hg38.exon.sort.bedtools.merge.bed > hg38.intron.bed
# 3.
cat hg38.intron.bed  | sort -k1,1V -k2,2n >hg38.intron.sort.bed
# 4.
bedtools merge -i hg38.intron.sort.bed > hg38.intron.sort.bedtools.merge.bed
# 5.
rm hg38.gene.sort.bed
rm hg38.intron.bed
rm hg38.intron.sort.bed
# 1.
zcat gencode.v44.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sort -k1,1V -k2,2n >hg38.gene.sort.bed
# 2.size
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg38.chromInfo"  > hg38.genome
# 3.size
sed '1d' hg38.genome > hg38_1.genome
# 4.
cat hg38_1.genome|sort -k1,1V -k2,2n >hg38.sort.genome
5.bedtools complement
bedtools complement -i hg38.gene.sort.bed -g hg38.sort.genome > hg38.intergenetic.sort.bedtools.bed
# 6.
rm hg38.gene.sort.bed
rm hg38.genome
rm hg38_1.genome
# UCSC,Tss0-2000bp，，+TSS2000bp，-Tss2000bp