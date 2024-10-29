import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve, PrecisionRecallDisplay
from sklearn import metrics
import os.path

"""
All helper functions for running the model and plotting the results.
"""
def abs_log1p(x):
    return np.log2(abs(x) + 1)


def chrom_train_test_split(pos, neg, chroms, high_pips, low_pips, class_weight=None):
    """
    Leave one chromosome out training and testing on high/low data only.
    """
    y_real_all = []
    y_prob_all = []
    peak_gene_order = []
    np.random.seed(1)
    for chr_str in chroms:
        print(chr_str)
        train_pos = high_pips.chr != chr_str
        train_neg = low_pips.chr != chr_str
        X_train = pd.concat([pos[train_pos], neg[train_neg]])
        X_test = pd.concat([pos[~train_pos], neg[~train_neg]])
        peak_gene_order = peak_gene_order + list(high_pips[~train_pos].index.values)
        peak_gene_order = peak_gene_order + list(low_pips[~train_neg].index.values)

        y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
        y_test = np.array([1] * sum(train_pos == 0) + [0] * sum(train_neg == 0))
        regr = LogisticRegression(
            random_state=1, max_iter=1000, class_weight=class_weight
        )

        # train the predictor
        regr.fit(X_train, y_train)
        y_prob = regr.predict_proba(X_test)
        y_prob = y_prob[:, 1]

        # save the prob
        y_real_all = y_real_all + list(y_test)
        y_prob_all = y_prob_all + list(y_prob)

    return y_real_all, y_prob_all, peak_gene_order


# Full model
def run_full_lr(high_pips, low_pips, annotation_cols, var_annots_to_pred):
    """
    Leave one chr out training, but predict on full data.
    """
    high_pip_pos = high_pips.loc[:, annotation_cols].fillna(0)
    low_pips_neg = low_pips.loc[:, annotation_cols].fillna(0)
    # and turn the log1p
    pos = high_pip_pos.apply(lambda x: abs_log1p(x))
    neg = low_pips_neg.apply(lambda x: abs_log1p(x))

    prob_dfs = {}

    np.random.seed(1)
    for chr_str in var_annots_to_pred.chr.unique():
        print(chr_str)
        train_pos = high_pips.chr != chr_str
        train_neg = low_pips.chr != chr_str
        X_train = pd.concat([pos[train_pos], neg[train_neg]])

        y_train = np.array([1] * sum(train_pos == 1) + [0] * sum(train_neg == 1))
        regr = LogisticRegression(
            random_state=1, max_iter=1000, class_weight="balanced"
        )

        # train the predictor
        regr.fit(X_train, y_train)

        # predict on all chr variants/peaks
        chr_data = var_annots_to_pred[var_annots_to_pred.chr == chr_str]
        chr_annots = (chr_data.loc[:, annotation_cols].fillna(0).apply(lambda x: abs_log1p(x)))

        # save the betas from chrom20
        if chr_str == "chr20":
            betas_df = pd.DataFrame(
                data={"betas": regr.coef_[0], "feature": chr_annots.columns}
            )

        prob_df = pd.DataFrame(
            {
                "peak_name": chr_data.peak_name,
                "phenotype_id": chr_data.phenotype_id,
                "predict_prob_high": regr.predict_proba(chr_annots)[:, 1],
            }
        )
        prob_dfs[chr_str] = prob_df

        all_peak_probs = pd.concat(prob_dfs)
        all_annots_with_predictions = pd.concat([var_annots_to_pred.set_index(["peak_name", "phenotype_id"]),
                                                 all_peak_probs.set_index(["peak_name", "phenotype_id"])],
                                                 axis=1)

    try:
        betas_df
    except NameError:
        betas_df = pd.DataFrame(
                data={"betas": regr.coef_[0], "feature": chr_annots.columns}
            )

    return all_annots_with_predictions, betas_df, prob_dfs, regr


def make_enrichment_plot(annots_with_pred_bin, annotations, numeric_columns=None):
    bins = [0, 0.01, 0.1, 0.5, 0.9, 1]
    labels = ["PIP<0.01", "0.01<PIP<0.1", "0.1<PIP<0.5", "0.5<PIP<0.9", "0.9<PIP"]
    annots_with_pred_bin["pip_bin"] = pd.cut(
        annots_with_pred_bin.pip,
        bins=bins,
        labels=labels,
        right=True,
        include_lowest=True,
    )

    mean_arr = pd.DataFrame(0.0, index=annotations, columns=labels)
    tot_arr = pd.DataFrame(0.0, index=annotations, columns=labels)

    for i, conseq in enumerate(annotations):
        for j, label in enumerate(labels):
            if conseq == "mean_start_distance":
                mean_arr.at[conseq, label] = (
                    annots_with_pred_bin.query("pip_bin==@label")
                    .groupby(["phenotype_id"], as_index=True)[conseq]
                    .agg("mean")
                    .mean()
                )
                continue
            elif numeric_columns is not None and conseq in numeric_columns:
                mean_arr.at[conseq, label] = (
                    annots_with_pred_bin.query("pip_bin==@label")
                    .groupby(["phenotype_id"], as_index=True)[conseq]
                    .agg("max")
                    .mean()
                )
                continue
            annot_bin_df = (
                annots_with_pred_bin.query("pip_bin==@label")
                .groupby(["phenotype_id"], as_index=True)[annotations]
                .agg("any")
            )
            mean_arr.at[conseq, label] = annot_bin_df[
                conseq
            ].mean()  # ignores NaN by default
            tot_arr.at[conseq, label] = annot_bin_df[
                conseq
            ].sum()

            norm = mean_arr["PIP<0.01"].values
    FE = (mean_arr.T / norm).T
    annotations = FE.index.values
    FE_melt_new = FE.reset_index().melt(id_vars="index")
    tot_melt = tot_arr.reset_index().melt(id_vars='index')

    mean_arr.to_csv(f'mean_array_by_pip.tsv', sep='\t', header=True, index=True)
    tot_arr.to_csv(f'totals_array_by_pip.tsv', sep='\t', header=True, index=True)

    fig, ax = plt.subplots(figsize=(15, 6))
    palette = ["#929591", "#F97306", "#FE420F"]
    sns.set_palette(palette)
    sns.barplot(FE_melt_new, x="index", y="value", hue="variable")
    annotation_labels = [
        x.replace("_", " ")
        for x in annotations
    ]
    ax.set_xticklabels(annotation_labels, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Fold Enrichment", fontsize=25)
    ax.set_xlabel("")
    ax.legend(title="", fontsize=15)
    ax.set_title(f'Peak Enrichment by Max PIP, \n with Model Bin Groups', fontsize=25)
    for loc in np.arange(0, 16, 2.5):
        ax.axhline(loc, c="k", ls="--", lw=0.35, zorder=0)

    # clip axis
    ax.set_ylim(0, 8)
    fig.tight_layout()
    fig.savefig(f'peak_pip_enrichment_by_category.png', dpi=300)

    return mean_arr, FE, fig


def add_q_bins(all_annots_with_predictions):
    all_annots_with_predictions["q_bins"] = pd.qcut(all_annots_with_predictions.predict_prob_high, 4)
    all_annots_with_predictions["q_bins"] = (all_annots_with_predictions.q_bins.cat.rename_categories(
            [
                f"{i.left} to {i.right}"
                for i in all_annots_with_predictions["q_bins"].cat.categories
            ]
        )
    )
    all_annots_with_predictions["pip"] = all_annots_with_predictions.max_pip
    res = pd.concat([all_annots_with_predictions,
                    pd.get_dummies(all_annots_with_predictions.q_bins)],
                    axis=1)

    return res


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest="peak_gene_dfs", nargs='+', default=[])
    parser.add_argument("-a", dest="annotation_columns", nargs='+', default=[],
        help="binary annotation columns to use in predictions")
    parser.add_argument("--n", dest="numeric_columns", nargs='*', default=[],
        help="numeric annotation columns to use in predictions (treated differently)")
    parser.add_argument("-c", dest="column_peaks", type=str, required=True)
    parser.add_argument("-v",dest="vars_in_peaks",
        help="parquet with vars in peaks annotations",
        type=str, required=True,
    )
    args = parser.parse_args()

    vars_in_peaks = pd.read_parquet(args.vars_in_peaks)

    fm_list = []
    for file_path in args.peak_gene_dfs:
        fm_list.append(pd.read_parquet(file_path))

    peak_gene_df_filt = pd.concat(fm_list)
    if args.column_peaks in peak_gene_df_filt.columns:
        peak_gene_df_filt.drop(columns=[args.column_peaks], inplace=True)

    peak_gene_df_filt.to_parquet('combined_peak_gene_df.parquet')
    # NEW BOUNDARIES
    high_pips = peak_gene_df_filt.query("max_pip >= 0.5")
    low_pips = peak_gene_df_filt.query("max_pip <= 0.01")
    print("high pip shape: ", high_pips.shape, "low pip shape: ", low_pips.shape)
    print("filtered peak gene df shape: ", peak_gene_df_filt.shape)
    # reshape if necessary
    high_thresh = 0.4
    while (high_thresh >= 0.1) & (high_pips.shape[0] <= 10):
        print('Need more high pips for training data. Re-sizing')
        high_pips = peak_gene_df_filt.query(f"max_pip >= {high_thresh}")
        print(f'High threshold = {high_thresh}, size = {high_pips.shape}')
        high_thresh -= .1
    low_thresh = 0.02
    while (low_thresh <= 0.1) & (low_pips.shape[0] <= 10):
        print('Need more low pips for training data. Re-sizing')
        low_pips = peak_gene_df_filt.query(f"max_pip <= {low_thresh}")
        print(f'Low threshold = {low_thresh}, size = {low_pips.shape}')
        high_thresh += .01

    print("Now running the train test split to get accuracy")

    if args.numeric_columns is not None:
        annotation_cols = np.hstack(('mean_start_distance', args.annotation_columns, args.numeric_columns))
    else:
        annotation_cols = np.hstack(('mean_start_distance', args.annotation_columns))
    high_pip_pos = high_pips.loc[:, annotation_cols].fillna(0)
    low_pips_neg = low_pips.loc[:, annotation_cols].fillna(0)

    # and turn the log1p
    pos = high_pip_pos.apply(lambda x: abs_log1p(x))
    neg = low_pips_neg.apply(lambda x: abs_log1p(x))

    # run test with this data for accuracy
    y_real, y_prob, variant_order_model = chrom_train_test_split(
        pos,
        neg,
        peak_gene_df_filt.chr.unique(),
        high_pips,
        low_pips,
        class_weight="balanced",
    )
    precision, recall, thresholds = precision_recall_curve(y_real, y_prob)
    fig, ax = plt.subplots(figsize=(10,8))
    pr_display = PrecisionRecallDisplay(precision=precision, recall=recall).plot(ax=ax)
    ax.set_ylim(0)
    ax.set_title(f'Train/Test ROC Curve, High and Low Peaks', fontsize=30)
    fpr, tpr, thresholds = metrics.roc_curve(y_real, y_prob, pos_label=1)
    ax.text(x=.8,y=.8, s=f'AUC: {round(metrics.auc(fpr, tpr), 3)}', fontsize=20)
    fig.savefig(f'peak_predictor_roc.png', dpi=300)
    print("AUC: ", metrics.auc(fpr, tpr))

    # set up
    peak_gene_df_filt.reset_index(inplace=True)
    peak_gene_df_filt = peak_gene_df_filt.rename(
        columns={"level_0": "phenotype_id", "level_1": "peak_name"}
    )

    print("Running full model.")
    all_annots_with_predictions, betas_df, probs_dfs, regr = run_full_lr(
        high_pips, low_pips, annotation_cols, peak_gene_df_filt
    )

    all_annots_with_predictions = add_q_bins(all_annots_with_predictions)

    all_annots_with_predictions = all_annots_with_predictions.merge(
        vars_in_peaks.groupby(["peak_name", "phenotype_id"])
        .size()
        .rename("num_vars_in_peak"),
        left_index=True,
        right_index=True,
    )
    betas_df = pd.concat([pd.DataFrame([[regr.intercept_[0], 'intercept']], columns=betas_df.columns), betas_df], ignore_index=True)

    all_annots_with_predictions.reset_index().to_parquet(f'peak_preds.parquet')
    betas_df.to_parquet(f'peak_betas.parquet')

    print("Making enrichment plot")
    all_plot_annots = np.hstack((annotation_cols, all_annots_with_predictions.q_bins.cat.categories))
    mean_arr, FE, fig = make_enrichment_plot(all_annots_with_predictions, all_plot_annots, args.numeric_columns)

    print('Done.')

if __name__ == '__main__':
    main()