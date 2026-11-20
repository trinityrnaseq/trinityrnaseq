#!/usr/bin/env python3
"""Extract CDS sequences from a GENCODE *.pc_transcripts.fa.gz (header has CDS:start-end),
write >transcript_id;gene_id, collapsing identical CDS sequences within a gene."""
import gzip, sys

def records(path):
    h, s = None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if h is not None:
                    yield h, ''.join(s)
                h, s = line[1:], []
            else:
                s.append(line)
    if h is not None:
        yield h, ''.join(s)

def main(inp, outp):
    seen = set(); n = 0
    with open(outp, 'w') as out:
        for h, seq in records(inp):
            f = h.split('|')
            cds = [x for x in f if x.startswith('CDS:')]
            if not cds:
                continue
            a, b = map(int, cds[0][4:].split('-'))
            cseq = seq[a-1:b]
            key = (f[1], cseq)
            if key in seen:
                continue
            seen.add(key); n += 1
            out.write('>%s;%s\n' % (f[0], f[1]))
            for i in range(0, len(cseq), 60):
                out.write(cseq[i:i+60] + '\n')
    sys.stderr.write('%d unique CDS from %d genes\n' % (n, len({g for g, _ in seen})))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
