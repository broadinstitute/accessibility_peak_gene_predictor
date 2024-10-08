version 1.0

workflow peak_gene_predictor {
    input {
        File eQTL_annotations_file
        String peak_column
        Array[String] prediction_categories # dont need to include start distance in this, it's automatic.
        Array[String] remove_categories
        Array[String]? numeric_categories # if we have prediction categories that are not 0/1 but rather numeric (e.g. E2G #enhancers), include them here. Max value per peak-gene will be taken.
        Int distance_threshold = 250 # min distance allowed between var and start site
        String git_branch = "main"
    }

    call gather_variants_in_peaks_groups {
        input:
            eQTL_annotations_file=eQTL_annotations_file,
            peak_column=peak_column,
            distance_threshold = distance_threshold,
            git_branch=git_branch
    }

    scatter (groups in gather_variants_in_peaks_groups.groups_list) {
        call build_peak_gene_df {
            input:
                vars_in_peaks=gather_variants_in_peaks_groups.vars_in_peaks,
                groups = groups,
                prediction_categories=prediction_categories,
                remove_categories=remove_categories,
                numeric_categories=numeric_categories,
                git_branch=git_branch
        }
    }

    call gather_peak_genes_and_run_model {
        input:
            peak_gene_dfs = build_peak_gene_df.peak_gene_df,
            files_size = size(build_peak_gene_df.peak_gene_df, "GB"),
            prediction_categories=prediction_categories,
            numeric_categories=numeric_categories,
            vars_in_peaks=vars_in_peaks,
            peak_column = peak_column,
            git_branch=git_branch
    }

    output {
        File vars_in_peaks = gather_variants_in_peaks_groups.vars_in_peaks
        File peak_pip_plot = gather_peak_genes_and_run_model.peak_pip_plot
        File roc_curve = gather_peak_genes_and_run_model.roc_curve
        File model_beta_df = gather_peak_genes_and_run_model.model_beta_df
        File model_peak_prediction = gather_peak_genes_and_run_model.model_peak_prediction
        File mean_array_by_pip = gather_peak_genes_and_run_model.mean_array_by_pip
        File total_array_by_pip = gather_peak_genes_and_run_model.total_array_by_pip
    }
}

task gather_variants_in_peaks_groups {
    input {
        File eQTL_annotations_file
        String peak_column
        Int distance_threshold
        String git_branch = "main"
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/accessibility_peak_gene_predictor.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/peak_gene_predictor/vars_in_peaks_pairs.py -f ${eQTL_annotations_file} -p ${peak_column} -d ${distance_threshold}
    }

    output {
        File vars_in_peaks = "vars_in_peaks.parquet"
        Array[File] groups_list = glob("group_file_*.txt")
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 2
        memory: "16GB"
        preemptible: 1
    }
}

task build_peak_gene_df {
    input {
        File vars_in_peaks
        File groups
        Array[String] prediction_categories
        Array[String] remove_categories
        Array[String]? numeric_categories
        String git_branch = "main"
    }
    String numeric_pre = if defined(numeric_categories) then "--n " else ""

    command {
        set -ex
        (git clone https://github.com/broadinstitute/accessibility_peak_gene_predictor.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/peak_gene_predictor/build_peak_gene_df.py -v ${vars_in_peaks} -g ${groups} -a ${sep=' ' prediction_categories} -r ${sep=' ' remove_categories} \
            ${numeric_pre}${sep=' ' numeric_categories}
    }

    output {
        File peak_gene_df = glob("*_peak_gene_df.parquet")[0]
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 4
        memory: "16GB"
    }
}

task gather_peak_genes_and_run_model {
    input {
        Array[String] peak_gene_dfs
        Float files_size
        Array[String] prediction_categories
        Array[String]? numeric_categories
        File vars_in_peaks
        String peak_column
        String git_branch = "main"
    }

    Int disk_size = ceil(files_size)
    String numeric_pre = if defined(numeric_categories) then "--n " else ""


    command {
        set -ex
        (git clone https://github.com/broadinstitute/accessibility_peak_gene_predictor.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/peak_gene_predictor/all_peak_gene_preds.py -p ${sep=' ' peak_gene_dfs} -a ${sep=' ' prediction_categories} -c ${peak_column} -v ${vars_in_peaks} \
            ${numeric_pre}${sep=' ' numeric_categories}
    }

    output {
        File peak_pip_plot = "peak_pip_enrichment_by_category.png"
        File roc_curve = "peak_predictor_roc.png"
        File model_beta_df = "peak_betas.parquet"
        File model_peak_prediction = "peak_preds.parquet"
        File mean_array_by_pip = "mean_array_by_pip.tsv"
        File total_array_by_pip = "totals_array_by_pip.tsv"
    }

    runtime {
        docker: 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        cpu: 4
        memory: "128GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}