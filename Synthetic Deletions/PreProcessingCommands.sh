echo "Starting preprocessing"

coverage=(5 10 20 30 40 50)
# coverage=(5 10)

for cov in ${coverage[@]}
do
    cd Raw_Seq_Files_CCS
    samtools faidx modRand.fasta
    pbsim --strategy wgs --method qshmm --qshmm ../../data/QSHMM-RSII.model --depth $cov --length-min 10000 --length-mean 15000 --length-max 20000 --accuracy-mean 0.99 --prefix modRand${cov}xCCS --genome modRand.fasta
    cd ..
    minimap2 -ax map-pb Raw_Seq_Files_CCS/refRand.fasta Raw_Seq_Files_CCS/modRand${cov}xCCS_0001.fastq > Alignment_Files_CCS/modRand${cov}xCCS.sam
    samtools view -@ 2 -Sb -o Alignment_Files_CCS/modRand${cov}xCCS.bam Alignment_Files_CCS/modRand${cov}xCCS.sam
    samtools sort -O bam -o Alignment_Files_CCS/sormodRand${cov}xCCS.bam Alignment_Files_CCS/modRand${cov}xCCS.bam
    samtools index Alignment_Files_CCS/sormodRand${cov}xCCS.bam
done

for cov in ${coverage[@]}
do
    cd Raw_Seq_Files_CCS
    mv modRand${cov}* ${cov}x
    cd ..
    cd Alignment_Files_CCS
    mv modRand${cov}* ${cov}x
    mv sormodRand${cov}* ${cov}x
    cd ..
done

echo "Preprocessing completed"

# cd Raw_Seq_Files_CCS
# samtools faidx modRand.fasta
# pbsim --strategy wgs --method qshmm --qshmm ../../data/QSHMM-RSII.model --depth 10 --length-min 10000 --length-mean 15000 --length-max 20000 --accuracy-mean 0.99 --prefix modRand10xCCS --genome modRand.fasta
# cd ..
# mkdir Alignment_Files_CCS
# minimap2 -ax map-pb Raw_Seq_Files_CCS/refRand.fasta Raw_Seq_Files_CCS/modRand10xCCS_0001.fastq > Alignment_Files_CCS/modRand10xCCS.sam
# samtools view -@ 2 -Sb -o Alignment_Files_CCS/modRand10xCCS.bam Alignment_Files_CCS/modRand10xCCS.sam
# samtools sort -O bam -o Alignment_Files_CCS/sormodRand10xCCS.bam Alignment_Files_CCS/modRand10xCCS.bam
# samtools index Alignment_Files_CCS/sormodRand10xCCS.bam

# cd Raw_Seq_Files_CCS
# mkdir 10x
# mv modRand10* 10x
# cd ..
# cd Alignment_Files_CCS
# mkdir 10x
# mv modRand10* 10x
# mv sormodRand10* 10x
# cd ..


















