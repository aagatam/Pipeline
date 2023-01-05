from Bio import SeqIO
import sys
#
input_fasta=sys.argv[1]
out_name=sys.argv[2]
# out_name='ES_test'
# input_fasta='//home/agata/Documents/Work/Dell_work/bisbee/Pipeline_test/Output/Test/genome/bisbee/exon_skip.altSeq.fasta'
output_fasta='.'.join([out_name,"final.fasta"])


seen = []
seen_ids = []
records = []

for record in SeqIO.parse(input_fasta, "fasta"):
    id = record.description.split(' ')[2].split('_')[:-1]
    id = '_'.join(id)
    if str(record.seq) not in seen:
        seen.append(str(record.seq))
        records.append(record)
        seen_ids.append(str(id))
    elif str(id) not in seen_ids:
        seen_ids.append(str(id))
        records.append(record)


#writing to a fasta file
SeqIO.write(records, output_fasta, "fasta")
