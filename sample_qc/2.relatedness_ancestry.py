import os
import hail as hl
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
from typing import List, Tuple
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("unified_sample_qc_a")
logger.setLevel(logging.INFO)

project_root = Path(__file__).parent.parent
print(project_root)

s3credentials = os.path.join(
    project_root, "hail_configuration_files/s3_credentials.json")
print(s3credentials)

storage = os.path.join(project_root, "hail_configuration_files/storage.json")

thresholds = os.path.join(
    project_root, "hail_configuration_files/thresholds.json")

with open(f"{s3credentials}", 'r') as f:
    credentials = json.load(f)

with open(f"{storage}", 'r') as f:
    storage = json.load(f)

with open(f"{thresholds}", 'r') as f:
    thresholds = json.load(f)


def get_related_samples_to_drop(rank_table: hl.Table, relatedness_ht: hl.Table) -> hl.Table:
    """
    Use the maximal independence function in Hail to intelligently prune clusters of related individuals, removing
    less desirable samples while maximizing the number of unrelated individuals kept in the sample set

    :param Table rank_table: Table with ranking annotations across exomes and genomes, computed via make_rank_file()
    :param Table relatedness_ht: Table with kinship coefficient annotations computed via pc_relate()
    :return: Table containing sample IDs ('s') to be pruned from the combined exome and genome sample set
    :rtype: Table
    """
    # Define maximal independent set, using rank list
    related_pairs = relatedness_ht.filter(
        relatedness_ht.kin > 0.08838835).select('i', 'j')
    n_related_samples = hl.eval(hl.len(
        related_pairs.aggregate(
            hl.agg.explode(
                lambda x: hl.agg.collect_as_set(x),
                [related_pairs.i, related_pairs.j]
            ),
            _localize=False)
    ))
    logger.info(
        '{} samples with at least 2nd-degree relatedness found in callset'.format(n_related_samples))
    max_rank = rank_table.count()
    related_pairs = related_pairs.annotate(id1_rank=hl.struct(id=related_pairs.i, rank=rank_table[related_pairs.i].rank),
                                           id2_rank=hl.struct(
                                               id=related_pairs.j, rank=rank_table[related_pairs.j].rank)
                                           ).select('id1_rank', 'id2_rank')

    def tie_breaker(l, r):
        return hl.or_else(l.rank, max_rank + 1) - hl.or_else(r.rank, max_rank + 1)

    related_samples_to_drop_ranked = hl.maximal_independent_set(related_pairs.id1_rank, related_pairs.id2_rank,
                                                                keep=False, tie_breaker=tie_breaker)
    return related_samples_to_drop_ranked.select(**related_samples_to_drop_ranked.node.id).key_by('data_type', 's')


def split_mt_by_relatedness(pruned_mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    """
    Split gnomAD Mt into unrelated and related MTs (based on hard-coded data paths)
    """
    related_samples_to_drop_ranked = hl.read_table(
        qc_temp_data_prefix('joint') + '.related_samples_to_drop.ht')
    rank_table = hl.import_table(rank_annotations_path(
        'joint'), impute=True).key_by('data_type', 's')
    pruned_mt = pruned_mt.annotate_cols(
        drop_unrelated_ranked=hl.is_defined(
            related_samples_to_drop_ranked[pruned_mt.col_key]),
        **rank_table[pruned_mt.col_key])
    pca_mt = pruned_mt.filter_cols(
        pruned_mt.drop_unrelated_ranked, keep=False).annotate_cols(related=False)
    related_mt = pruned_mt.filter_cols(
        pruned_mt.drop_unrelated_ranked, keep=True).annotate_cols(related=True)
    return pca_mt, related_mt


if __name__ == "__main__":
    # need to create spark cluster first before intiialising hail
    sc = pyspark.SparkContext()
    # Define the hail persistent storage directory
    tmp_dir = "hdfs://spark-master:9820/"
    temp_dir = os.path.join(os.environ["HAIL_HOME"], "tmp")
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")
    # s3 credentials required for user to access the datasets in farm flexible compute s3 environment
    # you may use your own here from your .s3fg file in your home directory
    hadoop_config = sc._jsc.hadoopConfiguration()

    hadoop_config.set("fs.s3a.access.key", credentials["mer"]["access_key"])
    hadoop_config.set("fs.s3a.secret.key", credentials["mer"]["secret_key"])

    #####################################################################
    ###################### INPUT DATA  ##############################
    #####################################################################
    CHROMOSOME = "WGS"
    # mt = hl.read_matrix_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_sex_annotations.mt")

    # ld pruning
    #pruned_ht = hl.ld_prune(mt.GT, r2=0.1)
    #pruned_mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))
    # pruned_mt.write(
    #    f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt", overwrite=True)
    pruned_mt = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt")
    # PC relate
    pruned_mt = pruned_mt.select_entries(
        GT=hl.unphased_diploid_gt_index_call(mt1.GT.n_alt_alleles()))

    eig, scores, _ = hl.hwe_normalized_pca(
        pruned_mt.GT, k=10, compute_loadings=False)
    scores.write(
        f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_pruned.pca_scores.ht", overwrite=True)

    relatedness_ht = hl.pc_relate(pruned_mt.GT, min_individual_maf=0.05,
                                  scores_expr=scores[pruned_mt.col_key].scores, block_size=4096, min_kinship=0.05, statistics='kin2')
    relatedness_ht.write(
        f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_relatedness.ht", overwrite=True)
    pairs = relatedness_ht.filter(relatedness_ht['kin'] > 0.125)
    related_samples_to_remove = hl.maximal_independent_set(
        pairs.i, pairs.j, keep=False)
    related_samples_to_remove.write(
        f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_related_samples_to_remove.ht", overwrite=True)

    pca_mt = pruned_mt.filter_cols(hl.is_defined(
        related_samples_to_remove[pruned_mt.col_key]), keep=False)
    related_mt = pruned_mt.filter_cols(hl.is_defined(
        related_samples_to_remove[pruned_mt.col_key]), keep=True)

    variants, samples = pca_mt.count()

    print(f"{samples} samples after relatedness step.")

    # Population pca

    plink_mt = pca_mt.annotate_cols(
        uid=pca_mt.s).key_cols_by('uid')
    hl.export_plink(plink_mt, f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_unrelated.plink",
                    fam_id=plink_mt.uid, ind_id=plink_mt.uid)
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(
        pca_mt.GT, k=20, compute_loadings=True)
    pca_af_ht = pca_mt.annotate_rows(
        pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(
        pca_af=pca_af_ht[pca_loadings.key].pca_af)
    pca_scores.write(
        f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_pca_scores.ht", overwrite=True)
    pca_loadings.write(
        f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_pca_loadings.ht", overwrite=True)

    pca_mt = pca_mt.annotate_cols(scores=pca_scores[pca_mt.col_key].scores)

    variants, samples = related_mt.count()
    print('Projecting population PCs for {} related samples...'.format(samples))
    related_scores = pc_project(related_mt, pca_loadings)
    relateds = related_mt.cols()
    relateds = relateds.annotate(scores=related_scores[relateds.key].scores)

    pca_mt.write(
        f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_pca.mt", overwrite=True)
    p = hl.plot.scatter(pca_mt.scores[0],
                        pca_mt.scores[1],
                        title='PCA', xlabel='PC1', ylabel='PC2')
    output_file(f"{temp_dir}/ddd-elgh-ukbb/pca.html")
    save(p)
