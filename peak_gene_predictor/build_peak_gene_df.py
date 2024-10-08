import pandas as pd
import numpy as np
import argparse
from ast import literal_eval as make_tuple


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v",dest="vars_in_peaks",
        help="parquet with vars in peaks annotations",
        type=str, required=True,
    )
    parser.add_argument("-g", dest='group_file',
        help="file containing peak-gene groups",
        type=str, required=True)
    parser.add_argument("-a", dest="annotation_columns", nargs='+', default=[])
    parser.add_argument("-r", dest="remove_columns", nargs='+', default=[])
    parser.add_argument("--n", dest="numeric_columns", nargs='*', default=[], required=False)
    args = parser.parse_args()

    group_file = open(args.group_file, 'r')
    Lines = group_file.readlines()
    groups = make_tuple(Lines[0])

    annotation_columns = args.annotation_columns
    remove_columns = args.remove_columns
    if args.numeric_columns is not None:
        numeric_columns = args.numeric_columns

    # first, we collect annotations in ALL of the 0/1 column types. numeric - these are calculated differently
    annotation_columns = np.hstack((annotation_columns, remove_columns))

    # read in the var data
    vars_in_peaks = pd.read_parquet(args.vars_in_peaks)
    vars_in_peaks.set_index(["phenotype_id", "peak_name"], inplace=True)
    vars_in_peaks.sort_index(inplace=True)

    # build a df
    # include numeric columns in the dataframe
    peak_gene_df = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(groups, names=('phenotype_id', 'peak_name')),
        columns=np.hstack(("max_pip", "mean_start_distance", annotation_columns, numeric_columns)) if args.numeric_columns is not None else np.hstack(("max_pip", "mean_start_distance", annotation_columns)),
    )
    for group in groups:
        temp_df = vars_in_peaks.loc[group]
        peak_gene_df.loc[group, annotation_columns] = (temp_df.loc[:, annotation_columns].any().astype(int))
        # max and mean for numeric cols should be the same if its E2G cols - others, take note that we take MAX value per group here.
        if args.numeric_columns is not None:
            peak_gene_df.loc[group, numeric_columns] = temp_df.loc[:, numeric_columns].max()
        peak_gene_df.loc[group, "mean_start_distance"] = temp_df.start_distance.mean()
        peak_gene_df.loc[group, "max_pip"] = temp_df.pip.max()
        peak_gene_df.loc[group, "chr"] = temp_df.chr.iloc[0]

    # now, remove any peaks that have an annotation in any of our remove categories
    peaks_to_keep = ~peak_gene_df[args.remove_columns].any(axis=1)
    peak_gene_df_filt = peak_gene_df.loc[peaks_to_keep]

    # print out which columns are resulting in most dropped
    print(peak_gene_df[args.remove_columns].sum(axis=1))

    fname = args.group_file
    name = fname.split('/')[-1].split('.')[0].split('_')[-1]
    peak_gene_df_filt.to_parquet(f'{name}_peak_gene_df.parquet')


if __name__ == '__main__':
    main()