# /usr/local/bin/trinityrnaseq/util/support_scripts/../../Trinity --single "/Users/foo/datasets/trinity/trinity_out_dir/read_partitions/Fb_0/CBin_0/c0.trinity.reads.fa" --output "/Users/foo/dataset
# s/trinity/trinity_out_dir/read_partitions/Fb_0/CBin_0/c0.trinity.reads.fa.out" --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --no_distributed_trinity_exec

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Parse Trinity commands')
    parser.add_argument('--commands', required=True)
    parser.add_argument('--output_dir', help='Output directory for rewriting input paths', required=True)
    parser.add_argument('--njobs', help='Number of jobs per node', required=True, type=int, default=10)
    args = parser.parse_args()
    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir += '/'
    input_commands = args.commands
    njobs = args.njobs
    inputs = []
    command_template = None
    with open(input_commands, 'rt') as f:
        with open('commands.txt', 'wt') as out:
            for line in f:
                line = line.strip()
                if line != '':
                    tokens = line.split(' ')
                    single_index = tokens.index('--single')
                    # rewrite input
                    input_file = tokens[single_index + 1]
                    input_file = output_dir + input_file[input_file.index("read_partitions"):]
                    input_file = input_file[:len(input_file) - 1]  # remove quote
                    inputs.append(input_file)
                    if command_template is None:
                        tokens[single_index + 1] = 'NA'
                        tokens[tokens.index('--output') + 1] = 'NA'
                        command_template = tokens

    with open("command_template.txt", "wt") as out:
        out.write(' '.join(command_template) + '\n')

    end = len(inputs)
    step = njobs
    shard = 0
    for i in range(0, end, step):
        offset_end = i + step
        offset_end = min(end, offset_end)
        # start index is inclusive, end index is exclusive
        with open("files-{}.txt".format(shard), "wt") as out:
            out.write('\n'.join(inputs[i:offset_end]))
        shard += 1
