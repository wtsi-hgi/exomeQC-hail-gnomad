# intersect 1000k genomes data with AKT to use for population PCA later on. 
#Results of overlaps are here: 
# /lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_AKT_intersections
DIR="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data"
akt="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/AKT_data/wes.hg38.vcf.gz"
outdir="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/intersections"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr10.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr11.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr12.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr13.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr14.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr15.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr16.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr17.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr18.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr19.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr1.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr20.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr21.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr22.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr2.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr3.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr4.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr5.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr6.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr7.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr8.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chr9.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chrX.recalibrated_variants.vcf.gz"
vcf="/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_data/CCDG_13607_B01_GRM_WGS_2019-02-19_chrY.recalibrated_variants.vcf.gz"

bcftools isec -p ${outdir}/chr10  -n~11 -c all -Oz ${vcf10} ${akt}
bcftools isec -p ${outdir}/chr10 -Oz ${vcf10} ${akt}

chr="chrY"

cp ${chr}/0000.vcf.gz.tbi ../../1kg_AKT_intersections/${chr}.vcf.gz.tbi
cp ${chr}/0000.vcf.gz ../../1kg_AKT_intersections/${chr}.vcf.gz
# final files:
/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/1kg_pca_hail/1kg_AKT_intersections
bcftools concat -o WES_AKT_1kg_intersection.vcf.bgz --output-type b chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz chrX.vcf.gz chrY.vcf.gz