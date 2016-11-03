from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint

human_reads = SeqIO.read(open("sequence.gb"), "genbank")

human_frags=[]
limit=len(human_reads.seq)
for i in range(0, 500):
    start = randint(0,limit-200)
    end = start+200
    human_frag = human_record.seq[start:end]
    record = SeqRecord(human_frag,'fragment_%i' % (i+1), '', '')
    human_frags.append(record)

output_handle = open("humanreads.fastq", "w")
SeqIO.write(human_frags, output_handle, "fastq")
output_handle.close()
