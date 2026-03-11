import pyfastx
import csv

K = 31

def rev_comp(seq):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]

def canonical(kmer):
    rc = rev_comp(kmer)
    return kmer if kmer < rc else rc

def get_kmers(seq, k=31):
    return [canonical(seq[i:i+k]) for i in range(len(seq) - k + 1)]

kmers_sig = set()
with open("/srv/home/mlef0011/VDARK/rawdata/kmer/tumour_specific.txt") as f:
    for line in f:
        kmer = line.strip().split()[0]
        kmers_sig.add(kmer)

kmers_canonical = {canonical(k) for k in kmers_sig}

output_file = "/srv/home/mlef0011/VDARK/rawdata/reads/read_kmer_association.tsv"

with open(output_file, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['read_ID', 'kmer'])
    
    for fq_file in [
        "/srv/home/mlef0011/VDARK/rawdata/reads/tumour_R1_tumour_specific.fq",
        "/srv/home/mlef0011/VDARK/rawdata/reads/tumour_R2_tumour_specific.fq"
    ]:
        for name, seq, qual in pyfastx.Fastq(fq_file, build_index=False):
            hits = {k for k in get_kmers(seq) if k in kmers_canonical}
            for kmer in hits:
                writer.writerow([name, kmer])

print("Done")