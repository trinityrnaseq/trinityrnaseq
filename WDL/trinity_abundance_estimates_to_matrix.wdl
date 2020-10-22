version 1.0

workflow trinity_abundance_estimates_to_matrix {
    input {

        Array[File] sample_results
        String output_base_name
        String input_format # RSEM|eXpress|kallisto|salmon
        String cross_sample_norm = "TMM" # TMM|UpperQuartile|none
        File? gene_to_transcript_map # gene-to-transcript mapping file
        Int preemptible = 2
        String memory = "2G"
        Int extra_disk_space = 4
        String docker = "trinityrnaseq/trinityrnaseq-wdl:2.11.0"
    }

    call abundance_estimates_to_matrix {
        input:
            sample_results = sample_results,
            output_base_name = output_base_name,
            input_format = input_format,
            cross_sample_norm = cross_sample_norm,
            gene_trans_map  = gene_to_transcript_map,
            extra_disk_space = extra_disk_space,
            preemptible = preemptible,
            memory = memory,
            docker = docker
    }

    output {
        File tpm = abundance_estimates_to_matrix.tpm
        File counts = abundance_estimates_to_matrix.counts
        Array[File] normalized_matrix = abundance_estimates_to_matrix.normalized_matrix
        Array[File] tmm_info = abundance_estimates_to_matrix.tmm_info
    }
}

task abundance_estimates_to_matrix {
    input {
        String input_format
        String cross_sample_norm
        Array[File] sample_results
        String output_base_name
        Int preemptible
        File? gene_trans_map
        String memory
        String docker
        Int extra_disk_space
    }
    String prefix = "~{output_base_name}"

    command <<<
        set -e

        gene_trans_map="~{gene_trans_map}"
        if [ "$gene_trans_map" == "" ]; then
            gene_trans_map="none"
        fi

        $TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
        --out_prefix ~{output_base_name} \
        --est ~{input_format} \
        --cross_sample_norm ~{cross_sample_norm} \
        --gene_trans_map $gene_trans_map \
        ~{sep=' ' sample_results}

        # remove .isoform from output names
        for f in *; do mv "$f" "${f/.isoform/}"; done;
    >>>

    output {
        File tpm = "~{prefix}.TPM.not_cross_norm"
        File counts = "~{prefix}.counts.matrix"
        Array[File] normalized_matrix = glob("~{prefix}.~{cross_sample_norm}.EXPR.matrix")
        Array[File] tmm_info = glob("~{prefix}.TPM.not_cross_norm.TMM_info.txt")
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(sample_results, "GB") + extra_disk_space) + " HDD"
        cpu: 1
        preemptible: preemptible
    }
}
