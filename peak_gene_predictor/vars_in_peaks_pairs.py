import pandas as pd
import numpy as np
import argparse
from itertools import islice

def chunks(data, SIZE=10000):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k for k in islice(it, SIZE)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f",dest="finemapped_annotations",
        help="parquet with all fm var annotations",
        type=str, required=True,
    )
    parser.add_argument("-p", dest="peak_column",
        help='name of column to define peaks',
        type=str, required=True)

    args = parser.parse_args()

    # read in fm annotation data
    fm_annots_parq = pd.read_parquet(args.finemapped_annotations)
    # should already be done, but just a failsafe
    fm_annots_parq.drop_duplicates(inplace=True)
    peak_column = args.peak_column

    # remove vars close to tss site
    var_annots_far = fm_annots_parq[fm_annots_parq.start_distance > 250]
    print(
        "number variants dropped due to start distance: ",
        fm_annots_parq.shape[0] - var_annots_far.shape[0],
    )
    # get variants in any DHS peak
    vars_in_peaks = var_annots_far.query(f"{peak_column} == 1")
    print(
        "number variants in an DHS peak (after distance removal): ",
        vars_in_peaks.shape[0],
    )

    print(
        "number of variants in DHS peak dropped due to start distance: ",
        fm_annots_parq.query(f"{peak_column} == 1").shape[0]
        - vars_in_peaks.shape[0],
    )

    vars_in_peaks.to_parquet(f'vars_in_peaks.parquet')
    peak_gene_groups = vars_in_peaks.groupby(["phenotype_id", "peak_name"]).groups.keys()

    # choose scatter size based on number of peak-gene pairs
    if len(peak_gene_groups) > 1e6:
        size = 1000
    elif len(peak_gene_groups) < 1000:
        size = 1
    else:
        size = 100

    with open('group_file.txt', 'w') as f:
        for item in chunks(peak_gene_groups, len(peak_gene_groups)//size):
            f.write(f'{item}\n')


if __name__ == '__main__':
    main()