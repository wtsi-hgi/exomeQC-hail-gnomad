import os
import hail as hl
import pandas as pd
import pyspark
import json
import sys
import re
from pathlib import Path
import logging
from typing import Any, Counter, List, Optional, Tuple, Union
from bokeh.plotting import output_file, save, show
from gnomad_ancestry import pc_project, run_pca_with_relateds

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


POP_NAMES = {
    "afr": "African/African-American",
    "ami": "Amish",
    "amr": "Latino",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "eur": "European",
    "fin": "Finnish",
    "mde": "Middle Eastern",
    "nfe": "Non-Finnish European",
    "oth": "Other",
    "sas": "South Asian",
    "uniform": "Uniform",
    "sas_non_consang": "South Asian (F < 0.05)",
    "consanguineous": "South Asian (F > 0.05)",
    "exac": "ExAC",
    "bgr": "Bulgarian (Eastern European)",
    "est": "Estonian",
    "gbr": "British",
    "nwe": "North-Western European",
    "seu": "Southern European",
    "swe": "Swedish",
    "kor": "Korean",
    "sgp": "Singaporean",
    "jpn": "Japanese",
    "oea": "Other East Asian",
    "oeu": "Other European",
    "onf": "Other Non-Finnish European",
    "unk": "Unknown",
}

POP_COLORS = {
    "afr": "#941494",
    "ami": "#FFC0CB",
    "amr": "#ED1E24",
    "asj": "#FF7F50",
    "eas": "#108C44",
    "eur": "#6AA5CD",
    "fin": "#002F6C",
    "mde": "#33CC33",
    "nfe": "#6AA5CD",
    "oth": "#ABB9B9",
    "sas": "#FF9912",
    "uniform": "pink",
    "consanguineous": "pink",
    "sas_non_consang": "orange",
    "exac": "gray",
    "bgr": "#66C2A5",
    "est": "black",
    "gbr": "#C60C30",
    "nwe": "#C60C30",
    "seu": "#3CA021",
    "swe": "purple",
    "kor": "#4891D9",
    "sgp": "darkred",
    "jpn": "#BC002D",
    "oea": "#108C44",
    "oeu": "#6AA5CD",
    "onf": "#6AA5CD",
    "unk": "#ABB9B9",
    "": "#ABB9B9",
}


def get_reference_genome(
    locus: Union[hl.expr.LocusExpression, hl.expr.IntervalExpression],
    add_sequence: bool = False,
) -> hl.ReferenceGenome:
    """
    Returns the reference genome associated with the input Locus expression
    :param locus: Input locus
    :param add_sequence: If set, the fasta sequence is added to the reference genome
    :return: Reference genome
    """
    if isinstance(locus, hl.expr.LocusExpression):
        ref = locus.dtype.reference_genome
    else:
        assert isinstance(locus, hl.expr.IntervalExpression)
        ref = locus.dtype.point_type.reference_genome
    if add_sequence:
        ref = add_reference_sequence(ref)
    return ref


def filter_to_autosomes(
    t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters the Table or MatrixTable to autosomes only.
    This assumes that the input contains a field named `locus` of type Locus
    :param t: Input MT/HT
    :return:  MT/HT autosomes
    """
    reference = get_reference_genome(t.locus)
    autosomes = hl.parse_locus_interval(
        f"{reference.contigs[0]}-{reference.contigs[21]}", reference_genome=reference
    )
    return hl.filter_intervals(t, [autosomes])


def assign_population_pcs(
    pop_pca_scores: Union[hl.Table, pd.DataFrame],
    pc_cols: Union[hl.expr.ArrayExpression, List[str]],
    known_col: str = "known_pop",
    fit: Any = None,  # Type should be RandomForestClassifier but we do not want to import sklearn.RandomForestClassifier outside
    seed: int = 42,
    prop_train: float = 0.8,
    n_estimators: int = 100,
    min_prob: float = 0.9,
    output_col: str = "pop",
    missing_label: str = "oth",
) -> Tuple[
    Union[hl.Table, pd.DataFrame], Any
]:  # 2nd element of the tuple should be RandomForestClassifier but we do not want to import sklearn.RandomForestClassifier outside
    """
    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.

    As input, this function can either take:

    - A Hail Table (typically the output of `hwe_normalized_pca`). In this case,
        - `pc_cols` should be an ArrayExpression of Floats where each element is one of the PCs to use.
        - A Hail Table will be returned as output
    - A Pandas DataFrame. In this case:
        - Each PC should be in a separate column and `pc_cols` is the list of all the columns containing the PCs to use.
        - A pandas DataFrame is returned as output

    .. note::

        If you have a Pandas Dataframe and have all PCs as an array in a single column, the
        `expand_pd_array_col` can be used to expand this column into multiple `PC` columns.

    :param pop_pc_pd: Input Hail Table or Pandas Dataframe
    :param pc_cols: Columns storing the PCs to use
    :param known_col: Column storing the known population labels
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param seed: Random seed
    :param prop_train: Proportion of known data used for training
    :param n_estimators: Number of trees to use in the RF model
    :param min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param output_col: Output column storing the assigned population
    :param missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Hail Table or Pandas Dataframe (depending on input) containing sample IDs and imputed population labels, trained random forest model
    """
    from sklearn.ensemble import RandomForestClassifier

    hail_input = isinstance(pop_pca_scores, hl.Table)
    if hail_input:
        pop_pc_pd = pop_pca_scores.select(pca_scores=pc_cols).to_pandas()

        # Explode the PC array
        num_out_cols = min([len(x)
                            for x in pop_pc_pd["pca_scores"].values.tolist()])
        pc_cols = [f"PC{i+1}" for i in range(num_out_cols)]
        pop_pc_pd[pc_cols] = pd.DataFrame(pop_pc_pd["pca_scores"].values.tolist())[
            list(range(num_out_cols))
        ]

    else:
        pop_pc_pd = pop_pca_scores

    train_data = pop_pc_pd.loc[~pop_pc_pd[pca_scores].isnull()]

    N = len(train_data)

    # Split training data into subsamples for fitting and evaluating
    if not fit:
        random.seed(seed)
        train_subsample_ridx = random.sample(
            list(range(0, N)), int(N * prop_train))
        train_fit = train_data.iloc[train_subsample_ridx]
        fit_samples = [x for x in train_fit["s"]]
        evaluate_fit = train_data.loc[~train_data["s"].isin(fit_samples)]

        # Train RF
        #training_set_known_labels = train_fit[known_col].values
        training_set_known_labels = "unk"
        training_set_pcs = train_fit[pc_cols].values
        evaluation_set_pcs = evaluate_fit[pc_cols].values

        pop_clf = RandomForestClassifier(
            n_estimators=n_estimators, random_state=seed)
        pop_clf.fit(training_set_pcs, training_set_known_labels)
        print(
            "Random forest feature importances are as follows: {}".format(
                pop_clf.feature_importances_
            )
        )

        # Evaluate RF
        predictions = pop_clf.predict(evaluation_set_pcs)
        error_rate = 1 - sum(evaluate_fit[known_col] == predictions) / float(
            len(predictions)
        )
        print("Estimated error rate for RF model is {}".format(error_rate))
    else:
        pop_clf = fit

    # Classify data
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].values)
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].values)
    probs = pd.DataFrame(
        probs, columns=[f"prob_{p}" for p in pop_clf.classes_])
    pop_pc_pd = pd.concat([pop_pc_pd, probs], axis=1)
    probs["max"] = probs.max(axis=1)
    pop_pc_pd.loc[probs["max"] < min_prob, output_col] = missing_label
    pop_pc_pd = pop_pc_pd.drop(pc_cols, axis="columns")

    logger.info(
        "Found the following sample count after population assignment: {}".format(
            ", ".join(
                f"{pop}: {count}"
                for pop, count in Counter(pop_pc_pd[output_col]).items()
            )
        )
    )

    if hail_input:
        pops_ht = hl.Table.from_pandas(pop_pc_pd, key=list(pop_pca_scores.key))
        pops_ht.annotate_globals(
            assign_pops_from_pc_params=hl.struct(min_assignment_prob=min_prob)
        )
        return pops_ht, pop_clf
    else:
        return pop_pc_pd, pop_clf


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

    # mt = hl.read_matrix_table(
    #    f"{temp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_sex_annotations.mt")

    # ld pruning
    # pruned_ht = hl.ld_prune(mt.GT, r2=0.1)
    # pruned_mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))
    # pruned_mt.write(
    #    f"{tmp_dir}/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt", overwrite=True)
    pruned_mt = hl.read_matrix_table(
        f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ddd-elgh-ukbb/chr1_chr20_XY_ldpruned.mt")

    related_samples_to_drop = hl.read_table(
        f"{temp_dir}/ddd-elgh-ukbb/relatedness_ancestry/ddd-elgh-ukbb/chr1_chr20_XY_related_samples_to_remove.ht")

    logger.info("run_pca_with_relateds")
    pca_evals, pca_scores, pca_loadings = run_pca_with_relateds(
        pruned_mt, related_samples_to_drop)

    mt = pruned_mt.annotate_cols(scores=pca_scores[pruned_mt.col_key].scores)
    #mt = mt.annotate_cols(loadings=pca_loadings[pruned_mt.s].loadings)
    mt = mt.annotate_cols(known_pop="unk")
    #pca_scores = pca_scores.annotate(known_pop="unk")
    #pca_scores.write(f"{tmp_dir}/ddd-elgh-ukbb/pca_scores.ht", overwrite=True)
    #pca_loadings.write(f"{tmp_dir}/ddd-elgh-ukbb/pca_loadings.ht", overwrite=True)
    #pca_scores = hl.read_table(f"{temp_dir}/ddd-elgh-ukbb/pca_scores.ht")
    #pca_loadings = hl.read_table(f"{temp_dir}/ddd-elgh-ukbb/pca_loadings.ht")
    logger.info("assign population pcs")
   # population_assignment_table = assign_population_pcs(
    #    pca_scores, pca_loadings, known_col="known_pop")
    population_assignment_table = assign_population_pcs(
        mt.scores, pca_loadings, known_col="known_pop")
    population_assignment_table.write(
        f"{tmp_dir}/ddd-elgh-ukbb/pop_assignments.ht")
