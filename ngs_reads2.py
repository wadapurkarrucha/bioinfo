from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint

human_reads = SeqIO.read(open("sequence.gb"), "genbank")

human_frags=[]
limit=len(human_reads.seq)
for i in range(0, 100000):
    start = randint(0,limit-50)
    end = start+50
    human_frag = human_reads.seq[start:end]
    record = SeqRecord(human_frag,'fragment_%i' % (i+1), '', '')
    human_frags.append(record)

output_handle = open("humanreads.fastq", "w")
SeqIO.write(human_frags, output_handle, "fastq")
output_handle.close()
