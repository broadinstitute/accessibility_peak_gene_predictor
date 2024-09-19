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
    args = parser.parse_args()

    group_file = open(args.group_file, 'r')
    Lines = group_file.readlines()
    groups = make_tuple(Lines[0])
    annotation_columns = args.annotation_columns
    vars_in_peaks = pd.read_parquet(vars_in_peaks)
    vars_in_peaks.set_index(["phenotype_id", "peak_name"], inplace=True)
    vars_in_peaks.sort_index(inplace=True)

    # build a df
    peak_gene_df = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(groups, names=('phenotype_id', 'peak_name')),
        columns=np.hstack(("max_pip", "mean_start_distance", annotation_columns)),
    )
    for group in tqdm(groups):
        temp_df = vars_in_peaks.loc[group]
        peak_gene_df.loc[group, annotation_columns] = (temp_df.loc[:, annotation_columns].any().astype(int))
        peak_gene_df.loc[group, "mean_start_distance"] = temp_df.start_distance.mean()
        peak_gene_df.loc[group, "max_pip"] = temp_df.pip.max()
        peak_gene_df.loc[group, "chr"] = temp_df.chr[0]

    # just save it as the first name of the group pair in this set, it doesnt really matter
    name = '_'.join(list(groups)[0])
    peak_gene_df.to_parquet(f'{name}_peak_gene_df.parquet')


if __name__ == '__main__':
    main()