version 1.0

workflow trinity_differential_expression {
    input {
        File samples
        File counts
        String output_base_name
        String method = "DESeq2"
        Int preemptible = 2
        String memory = "2G"
        Int extra_disk_space = 4
        String docker = "trinityrnaseq/trinityrnaseq-wdl:2.11.0"
    }

    call differential_expression {
        input:
            matrix = counts,
            output_base_name = output_base_name,
            samples  = samples,
            method = method,
            extra_disk_space = extra_disk_space,
            preemptible = preemptible,
            memory = memory,
            docker = docker
    }

    output {
        File differential_expression_results = differential_expression.results
        File differential_expression_plot = differential_expression.pdf
    }
}

task differential_expression {
    input {
        File matrix
        Int extra_disk_space
        String output_base_name
        File samples
        String method
        Int preemptible
        String memory
        String docker
    }
    command <<<
        set -e

        $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
        -o . \
        --method ~{method} \
        --matrix ~{matrix} \
        --samples_file ~{samples}
    >>>

    output {
        File results = glob("*.DE_results")[0]
        File pdf = glob("*.pdf")[0]
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(matrix, "GB") + extra_disk_space) + " HDD"
        cpu: 1
        preemptible: preemptible
    }
}
