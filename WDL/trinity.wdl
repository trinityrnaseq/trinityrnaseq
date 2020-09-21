version development

workflow trinity {
    input {
        Array[File] left
        Array[File] right

        Int cpu_read_clustering_phase = 6
        String memory_read_clustering_phase = "150G"
        Int extra_disk_space_read_clustering_phase = 90
        Float disk_space_multiplier_read_clustering_phase = 6
        Int preemptible_read_clustering_phase = 2
        String extra_args = "--seqType fq --no_distributed_trinity_exec"

        Int jobs_per_node_assembly_phase = 40
        Int preemptible_assembly_phase = 2
        Int cpu_assembly_phase = 2
        String memory_assembly_phase = "5G"

        String memory_gather_fasta_phase = "1G"
        Int preemptible_gather_fasta_phase = 2

        Int preemptible = 2
        String docker = "trinityrnaseq/trinityrnaseq-wdl:1.0.0"
    }

    String run_id = "trinity_out_dir"

    call trinity_read_clustering {
        input:
            cpu = cpu_read_clustering_phase,
            memory = memory_read_clustering_phase,
            extra_disk_space = extra_disk_space_read_clustering_phase,
            disk_space_multiplier = disk_space_multiplier_read_clustering_phase,
            preemptible = preemptible_read_clustering_phase,
            docker = docker,
            run_id = run_id,
            extra_args = extra_args,
            left = left,
            right = right
    }

    call parse_read_clustering_commands {
        input:
            preemptible = preemptible,
            docker = docker,
            commands = trinity_read_clustering.commands,
            output_directory = trinity_read_clustering.output_directory, # _ + '/' + run_id + '/'
            jobs = jobs_per_node_assembly_phase
    }

    scatter(fasta_shard in parse_read_clustering_commands.fasta_shards) {
        call trinity_assemble {
            input:
                preemptible = preemptible_assembly_phase,
                docker = docker,
                cpu = cpu_assembly_phase,
                memory = memory_assembly_phase,
                input_files = read_lines(fasta_shard),
                command_template = parse_read_clustering_commands.command_template
        }
    }

    call gather_fastas {
        input:
            fastas=trinity_assemble.fasta_shard,
            preemptible = preemptible_gather_fasta_phase,
            docker = docker,
            cpu = 1,
            memory = memory_gather_fasta_phase
    }

    output {
        File fasta = gather_fastas.fasta
        File gene_trans_map = gather_fastas.gene_trans_map
    }

}

task gather_fastas {
    input {
        String docker
        Int cpu
        String memory
        Int preemptible
        Array[File] fastas
    }
    Int disk_space = ceil(1+size(fastas, "GB")*3)

    command <<<
        set -e

#        /software/monitor_script.sh &

        output_name="Trinity.fasta"
        input_files="~{sep="," fastas}"
        IFS=','
        read -ra files <<< "$input_files"
        nfiles=${#files[@]}
        for (( i=0; i<${nfiles}; i++ )); do
            cat ${files[$i]} >> $output_name
        done

        /usr/local/bin/trinityrnaseq/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map
    >>>

    output {
        File gene_trans_map = "Trinity.fasta.gene_trans_map"
        File fasta = "Trinity.fasta"
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

task trinity_assemble {
    input {
        String docker
        Int cpu
        String memory
        Int preemptible
        File command_template
        Array[File] input_files
    }
    Int disk_space = ceil(size(input_files, "GB")*3)

    command <<<
        set -e

#        /software/monitor_script.sh &

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

task parse_read_clustering_commands {
    input {
        String docker
        String output_directory
        File commands
        Int preemptible
        Int jobs
    }

    command <<<
        set -e

        python /software/parse_commands.py --commands ~{commands} --output_dir ~{output_directory} --njobs ~{jobs}
    >>>

    output {
        Array[File] fasta_shards = glob("files-*.txt")
        File command_template = "command_template.txt"
    }

    runtime {
        docker: "~{docker}"
        memory: "1GB"
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(commands, "GB")*3) + " HDD"
        cpu: 1
        preemptible: preemptible
    }
}

task trinity_read_clustering {
    input {
        String docker
        Int cpu
        String memory
        String run_id
        Int extra_disk_space
        Float disk_space_multiplier
        Int preemptible
        Array[File] left
        Array[File] right
        String extra_args
    }
    Int disk_space = ceil((size(left, "GB")+size(right, "GB"))*disk_space_multiplier + extra_disk_space)
    command <<<
        set -e

#        /software/monitor_script.sh &

        mkdir ~{run_id}

        command_mem=$(cat /proc/meminfo | grep MemAvailable | awk 'BEGIN { FS=" " } ; { print int($2) }')
        command_mem=$(($command_mem/1000000))
        command_mem=$command_mem"G"

        Trinity --left ~{sep=',' left} --right ~{sep=',' right} --max_memory ${command_mem} --CPU $(nproc) --output `pwd`/~{run_id} ~{extra_args}

    >>>

    output {
        Directory output_directory = "~{run_id}"
        File commands = "~{run_id}/recursive_trinity.cmds"
    }

    runtime {
        docker: "~{docker}"
        memory: "~{memory}"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

