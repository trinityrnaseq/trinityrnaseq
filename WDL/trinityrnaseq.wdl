version 1.0

workflow trinityrnaseq {
    input {
        Array[File]? left
        Array[File]? right
        String output_base_name

        File? genome_guided_bam
        Int? genome_guided_max_intron

        Int cpu_read_clustering_phase = 6
        String memory_read_clustering_phase = "150G"
        Int disk_space_read_clustering_phase = 200
        Int preemptible_read_clustering_phase = 2
        String? seq_type
        String? extra_args # e.g. --min_contig_length 300 --jaccard_clip

        Int jobs_per_node_assembly_phase = 40
        Int preemptible_assembly_phase = 2
        Int cpu_assembly_phase = 2
        Float disk_space_multiplier_assembly_phase = 3
        String memory_assembly_phase = "5G"

        String memory_gather_fasta_phase = "1G"
        Int preemptible_gather_fasta_phase = 2

        Int cpu_genome_guided = 6
        String memory_genome_guided = "150G"
        Int extra_disk_space_genome_guided = 90
        Int disk_space_multiplier_genome_guided = 6
        Int preemptible_genome_guided = 2

        Int preemptible = 2
        String docker = "trinityrnaseq/trinityrnaseq-wdl:2.11.0"

    }
    String seq_type_use = select_first([seq_type, "fq"])
    String run_id = "trinity_out_dir"


    if(defined(genome_guided_bam)) {
        call trinity_genome_guided {
            input:
                run_id = run_id,
                extra_args = extra_args,
                cpu = cpu_genome_guided,
                memory = memory_genome_guided,
                extra_disk_space = extra_disk_space_genome_guided,
                disk_space_multiplier = disk_space_multiplier_genome_guided,
                preemptible = preemptible_genome_guided,
                genome_guided_bam = genome_guided_bam,
                genome_guided_max_intron = genome_guided_max_intron,
                output_base_name = output_base_name,
                docker = docker
        }
    }

    if(defined(left) || defined(right)) {
        call trinity_read_clustering {
            input:
                run_id = run_id,
                extra_args = extra_args,
                left = left,
                right = right,
                seq_type = seq_type_use,
                cpu = cpu_read_clustering_phase,
                memory = memory_read_clustering_phase,
                disk_space = disk_space_read_clustering_phase,
                preemptible = preemptible_read_clustering_phase,
                docker = docker
        }

        call create_shards {
            input:
                fastas = trinity_read_clustering.fastas,
                jobs = jobs_per_node_assembly_phase,
                preemptible = preemptible,
                docker = docker,
        }

        scatter(fasta_shard in create_shards.fasta_shards) {
            call trinity_assemble {
                input:
                    preemptible = preemptible_assembly_phase,
                    docker = docker,
                    cpu = cpu_assembly_phase,
                    disk_space_multiplier= disk_space_multiplier_assembly_phase,
                    memory = memory_assembly_phase,
                    input_files = read_lines(fasta_shard),
                    command_template = trinity_read_clustering.command_template
            }
        }

        call gather_fastas {
            input:
                fastas=trinity_assemble.fasta_shard,
                preemptible = preemptible_gather_fasta_phase,
                output_base_name = output_base_name,
                docker = docker,
                cpu = 1,
                memory = memory_gather_fasta_phase
        }
    }

    output {
        File? fasta = gather_fastas.fasta
        File? gene_trans_map = gather_fastas.gene_trans_map

        File? guided_fasta = trinity_genome_guided.fasta
        File? guided_gene_trans_map = trinity_genome_guided.gene_trans_map
    }

}

task gather_fastas {
    input {
        String docker
        Int cpu
        String memory
        Int preemptible
        String output_base_name
        Array[File] fastas
    }
    Int disk_space = ceil(1+size(fastas, "GB")*3)

    command <<<
        set -e

        output_name="~{output_base_name}.fasta"
        input_files="~{sep="," fastas}"
        IFS=','
        read -ra files <<< "$input_files"
        nfiles=${#files[@]}
        for (( i=0; i<${nfiles}; i++ )); do
            cat ${files[$i]} >> $output_name
        done

        /usr/local/bin/trinityrnaseq/util/support_scripts/get_Trinity_gene_to_trans_map.pl ~{output_base_name}.fasta > ~{output_base_name}.fasta.gene_trans_map
    >>>

    output {
        File gene_trans_map = "~{output_base_name}.fasta.gene_trans_map"
        File fasta = "~{output_base_name}.fasta"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

task trinity_assemble {
    input {
        String docker
        Int cpu
        String memory
        Int preemptible
        File command_template
        Array[File] input_files
        Float disk_space_multiplier
    }
    Int disk_space = ceil(size(input_files, "GB")*disk_space_multiplier)

    command <<<
        set -e

        python <<CODE

        import os

        input_files = "~{sep=',' input_files}".split(',')

        with open("~{command_template}", "rt") as f:
            command_template = f.readline().strip().split(' ')

        single_index = command_template.index('--single')
        output_index = command_template.index('--output')

        with open("commands.txt", "wt") as out:
            for i in range(len(input_files)):
                command_template[single_index + 1] = input_files[i]
                command_template[output_index + 1] = os.path.basename(input_files[i]) + ".out"
                out.write(" ".join(command_template) + "\n")
        CODE

        parallel --will-cite -a commands.txt --jobs $(nproc)

        find . -name "*out.Trinity.fasta" | /usr/local/bin/trinityrnaseq/util/support_scripts/partitioned_trinity_aggregator.pl --token_prefix TRINITY_DN --output_prefix Trinity
    >>>

    output {
        File fasta_shard = "Trinity.fasta"
    }

    runtime {
        docker: "~{docker}"
        memory: "~{memory}"
        bootDiskSizeGb: 12
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

task create_shards {
    input {
        String docker
        Array[String] fastas
        Int preemptible
        Int jobs
    }
    Int disk_space = 2

    command <<<
        set -e

        python <<CODE

        njobs = ~{jobs}
        fastas = "~{sep=',' fastas}".split(',')
        end = len(fastas)
        step = njobs
        shard = 0
        for i in range(0, end, step):
            offset_end = i + step
            offset_end = min(end, offset_end)
            # start index is inclusive, end index is exclusive
            with open("files-{}.txt".format(shard), "wt") as out:
                out.write('\n'.join(fastas[i:offset_end]))
            shard += 1
        CODE
    >>>

    output {
        Array[File] fasta_shards = glob("files-*.txt")
    }

    runtime {
        docker:docker
        memory: "1GB"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
    }
}

task trinity_genome_guided {
    input {
        String docker
        Int cpu
        String memory
        String run_id
        String output_base_name
        Int extra_disk_space
        Float disk_space_multiplier
        Int preemptible
        File? genome_guided_bam
        Int? genome_guided_max_intron
        String? extra_args
    }
    Int disk_space = ceil(disk_space_multiplier*size(genome_guided_bam, "GB") + extra_disk_space)

    command <<<
        set -e

        mkdir ~{run_id}

        command_mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2) }')
        command_mem=$(($command_mem/1000000))
        command_mem=$command_mem"G"

        Trinity \
        --max_memory $command_mem \
        --CPU $(nproc) \
        --output `pwd`/~{run_id} \
        ~{"--genome_guided_bam " + genome_guided_bam} \
        ~{"--genome_guided_max_intron " + genome_guided_max_intron} \
        ~{extra_args}

        mv ~{run_id}/Trinity-GG.fasta ~{output_base_name}.fasta
        mv ~{run_id}/Trinity-GG.fasta.gene_trans_map ~{output_base_name}.fasta.gene_trans_map

    >>>

    output {
        File fasta = "~{output_base_name}.fasta"
        File gene_trans_map = "~{output_base_name}.fasta.gene_trans_map"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

task trinity_read_clustering {
    input {
        String docker
        Int cpu
        String memory
        String run_id
        Int disk_space
        Int preemptible
        Array[File]? left
        Array[File]? right
        String? extra_args
        String seq_type
    }

#    Int disk_space = ceil((if defined(left) then disk_space_multiplier*size(left, "GB") else 0)+(if defined(right) then disk_space_multiplier*size(right, "GB") else 0) + extra_disk_space)
    String left_prefix = if(defined(left)) then "--left " else ""
    String right_prefix = if(defined(right)) then "--right " else ""

    command <<<
        set -e

        mkdir ~{run_id}

        command_mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2) }')
        command_mem=$(($command_mem/1000000))
        command_mem=$command_mem"G"

        Trinity \
        ~{left_prefix}~{sep=',' left} \
        ~{right_prefix}~{sep=',' right} \
        --seqType ~{seq_type} \
        --max_memory $command_mem \
        --CPU $(nproc) \
        --output `pwd`/~{run_id} \
        --no_distributed_trinity_exec \
        ~{extra_args}

        mkdir fasta
        python /software/parse_commands.py --commands ~{run_id}/recursive_trinity.cmds --output_dir fasta
    >>>

    output {
        Array[File] fastas = glob("fasta/*")
        File command_template = "command_template.txt"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

