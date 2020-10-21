# /usr/local/bin/trinityrnaseq/util/support_scripts/../../Trinity --single "/Users/foo/datasets/trinity/trinity_out_dir/read_partitions/Fb_0/CBin_0/c0.trinity.reads.fa" --output "/Users/foo/dataset
# s/trinity/trinity_out_dir/read_partitions/Fb_0/CBin_0/c0.trinity.reads.fa.out" --CPU 1 --max_memory 1G --run_as_paired --seqType fa --trinity_complete --full_cleanup --no_distributed_trinity_exec
import os
import shutil

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Parse Trinity commands')
    parser.add_argument('--commands', required=True)
    parser.add_argument('--output_dir', help='Output directory to move paritions to', required=True)
    args = parser.parse_args()
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    input_commands = args.commands
    inputs = []
    command_template = None

    with open(input_commands, 'rt') as f:
        with open('commands.txt', 'wt') as out:
            for line in f:
                line = line.strip()
                if line != '':
                    tokens = line.split(' ')
                    single_index = tokens.index('--single')
                    # get list of output files and command template
                    input_file = tokens[single_index + 1]
                    input_file = input_file.replace('"', '')
                    if not os.path.exists(input_file):
                        raise ValueError('{} not found'.format(input_file))
                    file_name = os.path.basename(input_file)
                    dest = os.path.join(output_dir, file_name)
                    if os.path.exists(dest):
                        base_file_name, extension = os.path.splitext(file_name)
                        counter = 1
                        new_name = '{}-{}{}'.format(base_file_name, counter, extension)
                        while os.path.exists(new_name):
                            counter += 1
                            new_name = '{}-{}{}'.format(base_file_name, counter, extension)
                        dest = os.path.join(output_dir, new_name)
                    shutil.move(input_file, dest)
                    inputs.append(dest)
                    if command_template is None:
                        tokens[single_index + 1] = 'NA'
                        tokens[tokens.index('--output') + 1] = 'NA'
                        command_template = tokens

    with open("command_template.txt", "wt") as out:
        out.write(' '.join(command_template) + '\n')
