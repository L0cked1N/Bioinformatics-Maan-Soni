# ori_finder.py
from Bio import SeqIO
from collections import defaultdict
import math

# ---------- Utilities ----------

def is_low_complexity(motif):
    return len(set(motif)) <= 2

def background_model(seq):
    total = len(seq)
    return {
        "A": seq.count("A") / total,
        "C": seq.count("C") / total,
        "G": seq.count("G") / total,
        "T": seq.count("T") / total
    }

# ---------- GC Skew ----------

def compute_gc_skew(seq):
    skew = [0]
    current = 0
    for base in seq:
        if base == "G":
            current += 1
        elif base == "C":
            current -= 1
        skew.append(current)
    return skew

def find_ori_by_gc_skew(seq):
    skew = compute_gc_skew(seq)
    return skew.index(min(skew))

# ---------- PWM (MEME core) ----------

def build_pwm(instances, k):
    pwm = []
    for i in range(k):
        col = {"A":1, "C":1, "G":1, "T":1}  # pseudocounts
        for motif in instances:
            col[motif[i]] += 1
        total = sum(col.values())
        for b in col:
            col[b] /= total
        pwm.append(col)
    return pwm

def score_with_pwm(seq, pwm, bg):
    score = 0
    k = len(pwm)
    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        llr = 0
        for j, base in enumerate(window):
            llr += math.log(pwm[j][base] / bg[base])
        score += llr
    return score

# ---------- MEME-style search ----------

def meme_like_search(seq, k):
    bg = background_model(seq)
    best_score = float("-inf")
    best_instances = None

    for i in range(len(seq) - k + 1):
        seed = seq[i:i+k]
        if is_low_complexity(seed):
            continue

        instances = [seed]

        for j in range(len(seq) - k + 1):
            window = seq[j:j+k]
            if not is_low_complexity(window):
                instances.append(window)

        pwm = build_pwm(instances, k)
        score = score_with_pwm(seq, pwm, bg)

        if score > best_score:
            best_score = score
            best_instances = instances

    return best_instances, best_score

def meme_k_selection(seq, k_min=6, k_max=15):
    results = []
    for k in range(k_min, k_max + 1):
        motifs, score = meme_like_search(seq, k)
        results.append((k, score, motifs))
    return max(results, key=lambda x: x[1])

# ---------- Final ORI Finder ----------

def find_ori_meme_style(fasta_file, region_size=500):
    record = SeqIO.read(fasta_file, "fasta")
    seq = str(record.seq).upper()

    ori_center = find_ori_by_gc_skew(seq)
    start = max(0, ori_center - region_size//2)
    end   = min(len(seq), ori_center + region_size//2)
    region = seq[start:end]

    k, score, motifs = meme_k_selection(region)

    return {
        "ori_position": ori_center,
        "ori_region": (start, end),
        "ori_sequence": region,
        "best_k": k,
        "log_likelihood": score,
        "motifs": motifs[:10]
    }
