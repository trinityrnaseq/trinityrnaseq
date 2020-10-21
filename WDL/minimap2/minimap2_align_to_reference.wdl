version 1.0

workflow minimap2_align_to_reference {
    input {
        File ref_fasta
        File assembled_fasta
        String output_base_name
        Int extra_disk_space = 1
        String flags = "asm20"
        String docker  = "trinityrnaseq/minimap2:2.17"
        Int cpu = 1
        String memory = "2G"
        Int preemptible = 2
    }

    call run_minimap2 {
        input:
            ref_fasta = ref_fasta,
            assembled_fasta = assembled_fasta,
            output_base_name = output_base_name,
            flags = flags,
            memory = memory,
            cpu = cpu,
            extra_disk_space = extra_disk_space,
            preemptible = preemptible,
            docker = docker
    }

    output {
        File bam = run_minimap2.bam
        File bam_index = run_minimap2.bam_index
    }
}

task run_minimap2 {
    input {
        File ref_fasta
        File assembled_fasta
        String output_base_name
        Int extra_disk_space
        String flags
        String docker
        Int cpu
        String memory
        Int preemptible
    }

    command <<<
        set -e

        minimap2 -a -cx ~{flags} --cs ~{ref_fasta} ~{assembled_fasta} > ~{output_base_name}.sam
        samtools sort ~{output_base_name}.sam > ~{output_base_name}.bam
        samtools index ~{output_base_name}.bam
    >>>

    output {
        File bam = "~{output_base_name}.bam"
        File bam_index = "~{output_base_name}.bam.bai"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(ref_fasta, "GB") + size(assembled_fasta, "GB") + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}
