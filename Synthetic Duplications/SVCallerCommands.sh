iteration=$1

echo "Starting SV calling"

coverage=(5 10 20 30 40 50)
# coverage=(5 10)

for cov in ${coverage[@]}
do
    echo "---------------Starting ${cov}x----------------"
    echo "Starting cuteSV for ${cov}x"
    mkdir modRand${cov}xCCSSVsCuteSV
    cuteSV ../Alignment_Files_CCS/${cov}x/sormodRand${cov}xCCS.bam ../Raw_Seq_Files_CCS/refRand.fasta modRand${cov}xCCSSVsCuteSV_$1.vcf modRand${cov}xCCSSVsCuteSV --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
    echo "CuteSV for ${cov}x completed"

    echo "-----------------------------------------------"
    
    echo "Starting Sniffles2 for ${cov}x"
    conda activate base
    sniffles -i ../Alignment_Files_CCS/${cov}x/sormodRand${cov}xCCS.bam -v modRand${cov}xCCSSVsSniffles_$1.vcf --reference ../Raw_Seq_Files_CCS/refRand.fasta
    echo "Sniffles2 for ${cov}x completed"

    echo "-----------------------------------------------"
    
    echo "Starting SVIM for ${cov}x"
    mkdir modRand${cov}xCCSSVsSVIM
    svim alignment modRand${cov}xCCSSVsSVIM ../Alignment_Files_CCS/${cov}x/sormodRand${cov}xCCS.bam ../Raw_Seq_Files_CCS/refRand.fasta
    mv modRand${cov}xCCSSVsSVIM/signatures/all.vcf .
    mv all.vcf modRand${cov}xCCSSVsSVIM_$1.vcf
    mv modRand${cov}x* ${cov}x
    echo "SVIM for ${cov}x completed"
    echo "---------------${cov}x completed----------------"
done

echo "SV calling completed"

echo "---------------Starting 5x----------------"
echo "Starting cuteSV for 5x"
mkdir modRand5xCCSSVsCuteSV
cuteSV ../Alignment_Files_CCS/5x/sormodRand5xCCS.bam ../Raw_Seq_Files_CCS/refRand.fasta modRand5xCCSSVsCuteSV.vcf modRand5xCCSSVsCuteSV --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
echo "CuteSV for 5x completed"

echo "-----------------------------------------------"

# echo "Starting Sniffles2 for 5x"
# conda activate base
# sniffles -i ../Alignment_Files_CCS/5x/sormodRand5xCCS.bam -v modRand5xCCSSVsSniffles.vcf --reference ../Raw_Seq_Files_CCS/refRand.fasta
# echo "Sniffles2 for 5x completed"

# echo "-----------------------------------------------"
    
# echo "Starting SVIM for 5x"
# mkdir modRand5xCCSSVsSVIM
# svim alignment modRand5xCCSSVsSVIM ../Alignment_Files_CCS/5x/sormodRand5xCCS.bam ../Raw_Seq_Files_CCS/refRand.fasta
# mv modRand5xCCSSVsSVIM/signatures/all.vcf .
# mv all.vcf modRand5xCCSSVsSVIM.vcf
# mkdir 5x
# mv modRand5x* 5x
# echo "SVIM for 5x completed"
# echo "---------------5x completed----------------"


# echo "---------------Starting 10x----------------"
# echo "Starting cuteSV for 10x"
# mkdir modRand10xCCSSVsCuteSV
# cuteSV ../Alignment_Files_CCS/10x/sormodRand10xCCS.bam ../Raw_Seq_Files_CCS/refRand.fasta modRand10xCCSSVsCuteSV.vcf modRand10xCCSSVsCuteSV --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5
# echo "CuteSV for 10x completed"

# echo "-----------------------------------------------"

# echo "Starting Sniffles2 for 10x"
# conda activate base
# sniffles -i ../Alignment_Files_CCS/10x/sormodRand10xCCS.bam -v modRand10xCCSSVsSniffles.vcf --reference ../Raw_Seq_Files_CCS/refRand.fasta
# echo "Sniffles2 for 10x completed"

# echo "-----------------------------------------------"
    
# echo "Starting SVIM for 10x"
# mkdir modRand10xCCSSVsSVIM
# svim alignment modRand10xCCSSVsSVIM ../Alignment_Files_CCS/10x/sormodRand10xCCS.bam ../Raw_Seq_Files_CCS/refRand.fasta
# mv modRand10xCCSSVsSVIM/signatures/all.vcf .
# mv all.vcf modRand10xCCSSVsSVIM.vcf
# mkdir 10x
# mv modRand10x* 10x
# echo "SVIM for 10x completed"
# echo "---------------10x completed----------------"




