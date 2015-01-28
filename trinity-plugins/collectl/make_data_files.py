__author__ = 'befulton'
from subprocess import call, Popen, PIPE
from collections import defaultdict
import os
import time
import re
import datetime
import sys

def total_seconds(td):
    # Since this function is not available in Python 2.6
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

def open_collectl(filename, start_time=None):
    bindir = os.path.dirname(sys.argv[0])
    args = [bindir+os.sep+"collectl", "-P", "-p", filename, "-sZ"]
    if start_time:
        midnight = datetime.datetime.combine(start_time, datetime.datetime.min.time())
        offsettime = total_seconds(midnight - start_time + datetime.timedelta(minutes=1))
        print "offsettime: %s s" % offsettime
        args += ["--offsettime", str(offsettime)]
    p = Popen(args, stdout=PIPE)
    return p

def wait_for_collectl(filename):
    finished = call(["gzip", "-tf", filename])
    count = 0
    while finished != 0 and count < 300:
        print "collectl has not yet finished"
        time.sleep(1)
        finished = call(["gzip", "-tf", filename])
        count += 1

def get_start_time(filename):
    p = open_collectl(filename)
    print p.stdout.readline()
    print p.stdout.readline()
    s = p.stdout.readline().split()

    p.stdout.close()
    return s[0:2]
    #return datetime.datetime.strptime(s[0] + ' ' + s[1], "%Y%m%d %H:%M:%S")

def replay(filename, start_time):
    if filename.endswith('gz'):
        p = open_collectl(filename, start_time)
        print p.stdout.readline()
        print p.stdout.readline()
        while True:
          retcode = p.poll() #returns None while subprocess is running
          line = p.stdout.readline()
          yield line
          if(retcode is not None):
            break
    else:
        with open(filename) as  f:
            for line in f:
                yield line

def timedelta(filename):
    print filename
    p = Popen(["zcat", filename], stdout=PIPE)
    p.stdout.readline()
    params = p.stdout.readline().split()
    i = next(param for param in params if param.startswith('-i'))
    p.stdout.close()
    return i[2:]


def prettyprocess(line):
    s = line.split()[29:]
    if not s:
        return None
    if s[0] in ("-bash", 'sh', '/bin/bash', 'bash', 'ln', '/bin/pwd', 'mkdir', 'date', 'touch', '/usr/bin/env'):  
        return None

    exes = ['fastool', 'ParaFly','Butterfly','ReadsToTranscripts', 'jellyfish',
             'inchworm', 'FastaToDeBruijn', 'QuantifyGraph', 'GraphFromFasta', 'CreateIwormFastaBundle',
             'bowtie-build', 'bowtie', 'Chrysalis', 'cat', 'sort', 'cp', 'wc', 'rm', 'find']
    perl_scripts = ['scaffold_iworm_contigs', 'Trinity', 'collectl', 'print_butterfly_assemblies',
              'partition_chrysalis_graphs_n_reads', 'fasta_filter_by_min_length', 'partitioned_trinity_aggregator']

    for k in exes:
        if s[0].endswith(k):
            return k

    if s[0] == 'samtools':
        return ('samtools_' + s[1]) if len(s) > 1 else 'samtools'

    if s[0] == '/bin/sort':
        return 'sort'

    if s[0] == 'java':
        if 'Butterfly.jar' in " ".join(s):
            return 'Butterfly'
        if 'ExitTester.jar' in " ".join(s):
            return 'ExitTester'
        if '-version' in " ".join(s):
            return 'java_version'

        return 'java'

    if s[0] == 'perl':
        for k in perl_scripts:
            if k in s[1]:
                return k
        return 'perl'

    if s[0] == '/usr/bin/perl' and 'collectl' in s[2]:
        return 'collectl'

    return os.path.basename(s[0])+'_unknown'


def build_datasets(line_generator):
    line_dict = defaultdict(list)
    sum_dict = defaultdict(lambda: [0]*27)
    grand_sum_dict = defaultdict(lambda: [0]*27)

    last_line = ''
    ordered_keys = []
    for line in line_generator:
        if line: last_line = line
        key = prettyprocess(line)
        if key:
            if key not in ordered_keys:
                ordered_keys.append(key)
            line_dict[key].append(line)
            s = line.split()
            data = sum_dict[key, s[0], s[1]]
            total = grand_sum_dict[tuple(s[0:2])]
            for i in range(27):
                try:
                    data[i] += float(s[i + 2])
                    total[i] += float(s[i + 2])
                except ValueError:
                    pass

    return ordered_keys, line_dict, sum_dict, grand_sum_dict, last_line.split()[0:2]

def write_times(start_time, end_time, colpar):
     ref = start_time[0] + " 00:00:00"
     e = datetime.datetime.strptime(" ".join(end_time), "%Y%m%d %H:%M:%S")
     r = datetime.datetime.strptime(ref, "%Y%m%d %H:%M:%S")

     with open("global.time", "w") as f:
        f.write("date %s\n" % " ".join(start_time))
        f.write("start %s\n" % ref)
        f.write("end %s\n" % " ".join(end_time))
        f.write("runtime %s\n" % total_seconds(e-r))
        f.write("interval %s\n" % colpar)


def write_files(keys, line_dict, sum_dict, grand_sum_dict):
    for id, tool in enumerate(keys):
        if line_dict[tool]:
            with open("%s.%s.data" % (id + 1, tool), "w") as f:
                f.writelines(line_dict[tool])
            with open("%s.%s.sum" % (id + 1, tool), "w") as f:
                lines = (key for key in sum_dict.keys() if key[0] == tool)
                for line in sorted(lines):
                    f.write(" ".join(list(line[1:]) + [('%f' % val).rstrip('0').rstrip('.')
                                                       for val in sum_dict[line]]) + "\n")
        with open("collectZ.proc", "w") as f:
            for k in sorted(grand_sum_dict.keys()):
                f.write(" ".join(list(k) + ['{0:g}'.format(val) for val in grand_sum_dict[k]]) + "\n")

if __name__ == '__main__':
    filename = next(file for file in os.listdir(".") if file.endswith(".gz"))
    wait_for_collectl(filename)
    colpar = timedelta(filename)
    start_time = get_start_time(filename)

    start = datetime.datetime.strptime(" ".join(start_time), "%Y%m%d %H:%M:%S")
    ordered_keys, line_dict, sum_dict, grand_sum_dict, end_time = build_datasets(replay(filename, start))
    write_files(ordered_keys, line_dict, sum_dict, grand_sum_dict)
    write_times(start_time, end_time, colpar)

    sum_files = sorted((file for file in os.listdir(".") if file.endswith(".sum")), key=lambda f: int(f.split('.')[0]))
    for id, sf in enumerate(sum_files):
        os.rename(sf, "%s.%s" % (id+1, sf))
