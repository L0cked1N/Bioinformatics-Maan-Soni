# ori_finder.py
from Bio import SeqIO
import statistics

def find_ori_single(fasta_file, window=200, step=50):
    record = SeqIO.read(fasta_file, "fasta")
    seq = str(record.seq).upper()

    skew_values = []
    positions = []

    for i in range(0, len(seq) - window, step):
        window_seq = seq[i:i+window]
        g = window_seq.count("G")
        c = window_seq.count("C")

        if g + c == 0:
            skew = 0
        else:
            skew = (g - c) / (g + c)

        skew_values.append(skew)
        positions.append(i)

    derivatives = []
    for i in range(len(skew_values) - 1):
        derivatives.append(abs(skew_values[i+1] - skew_values[i]))

    max_idx = derivatives.index(max(derivatives))
    ori_start = positions[max_idx]
    ori_end = ori_start + window

    return ori_start, ori_end


def find_ori_multi_scale(fasta_file, windows=[150,200,250,300], step=25):
    starts = []
    ends = []

    for w in windows:
        s, e = find_ori_single(fasta_file, window=w, step=step)
        starts.append(s)
        ends.append(e)

    # Consensus = median (more robust than mean)
    final_start = int(statistics.median(starts))
    final_end   = int(statistics.median(ends))

    record = SeqIO.read(fasta_file, "fasta")
    seq = str(record.seq).upper()
    ori_seq = seq[final_start:final_end]

    return ori_seq, final_start, final_end
