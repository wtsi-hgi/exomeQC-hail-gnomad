import logging
import random
from typing import Any, Counter, List, Optional, Tuple, Union

import hail as hl
import pandas as pd


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
        pop_pc_pd = pop_pca_scores.select(
            known_col, pca_scores=pc_cols).to_pandas()

        # Explode the PC array
        num_out_cols = min([len(x)
                            for x in pop_pc_pd["pca_scores"].values.tolist()])
        pc_cols = [f"PC{i+1}" for i in range(num_out_cols)]
        pop_pc_pd[pc_cols] = pd.DataFrame(pop_pc_pd["pca_scores"].values.tolist())[
            list(range(num_out_cols))
        ]

    else:
        pop_pc_pd = pop_pca_scores

    train_data = pop_pc_pd.loc[~pop_pc_pd[known_col].isnull()]

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
        training_set_known_labels = train_fit[known_col].values
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
