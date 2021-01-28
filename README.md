# sanger_gnomad_hail_qc

This repository contains hail scripts designed to import, analyse and execute sample_qc and variant_qc on a collection of data cohorts at the Sanger Institute.
The scripts are adapted from the gnomad sample and variant QC v2 and v3 analysis. (https://github.com/broadinstitute/gnomad_qc/blob/master/gnomad_qc/v2/ , https://github.com/broadinstitute/gnomad_qc/tree/master/gnomad_qc/v3).

Please visit https://confluence.sanger.ac.uk/display/HGI/WSI+Exome+Joint+Call for more information.
Please contact Pavlos (pa10) for any questions.

---

Python modules to install gnomad methods on hail cluster:
pip install cython
pip install gnomad

## Part 1: Import data to hail

### import_data/1.upload_vcf_from_farm_to_s3.sh

This script uploads data to s3 to be shared with hail cluster.

### import_data/2.import_vcf.py

This script imports the VCF files into hail and creates one hail matrixtable per chromosome.

## Part 2: Sample QC

### sample_qc_v3/1.hard_filters_sex_annotation.py

This script will apply gnomad hard filters and assign sex to each sample.

```python
usage: 1.hard_filters_sex_annotation.py [-h] [--matrixtable MATRIXTABLE]
                                        [--output-dir OUTPUT_DIR]

optional arguments:
  -h, --help            show this help message and exit

Input parameters:
  --matrixtable MATRIXTABLE
                        Full path of input matrixtable. chrX and chrY
                        variation should be included
  --output-dir OUTPUT_DIR
                        Full path of output folder to store results.
                        Preferably hdfs or secure lustre
```

### sample_qc_v3/2.ld_prune_relatednessPCA.py

### sample_qc_v3/3.population_PCA_prediction.py

### sample_qc_v3/4a.find_population_outliers_part1.py

### sample_qc_v3/4b.find_population_outliers_part2.py

### sample_qc_v3/5.filter_samples_fail_sample_qc.py

## Part 3: Variant QC

### variant_qc/1.generate_truthsets_final.py

### variant_qc/1a.family_stats_hail_de_novo.py

### variant_qc/2.create_rf_ht.py

### variant_qc/3.train_apply_finalise_RF.py

### variant_qc/4.rank_bin_plot_v2_denovo.py
