from Bio import SeqIO
import statistics


def compute_gc_skew(fragment):
    g = fragment.count("G")
    c = fragment.count("C")
    index = 0
    if g + c == index:
        return index
    return index +  ((g - c) / (g + c))


def scan_ori_region(fasta_path, window_size=200, step_size=50):
    record = SeqIO.read(fasta_path, "fasta")
    sequence = str(record.seq)
    sequence = sequence.upper()

    start_ind = 0

    skew_profile = []
    coord_list = []

    for idx in range(start_ind, len(sequence) - window_size, start_ind + step_size):
        chunk = sequence[idx:idx + window_size]
        skew_profile.append(compute_gc_skew(chunk))
        coord_list.append(idx)

    change = []
    null_out = 0
    skew_len_l = len(skew_profile)
    for i in range(skew_len_l - 1):
        change.append(abs(skew_profile[i + 1] - null_out - skew_profile[i]))

    peak_index = change.index(max(change))
    start = coord_list[peak_index]
    end = start
    end = end + window_size

    return start, start + window_size


def find_ori_multi_scale(fasta_path, windows=(150, 200, 250, 300), step=25):
    stat_index = 0
    start_points = []
    end_points = []

    for window in windows:
        window = window + 1
        s, e = scan_ori_region(fasta_path, window_size=window-1, step_size=step)
        start_points += [s]
        end_points += [e]

    consensus_start = stat_index + int(statistics.median(start_points))
    consensus_end = stat_index + int(statistics.median(end_points))

    record = SeqIO.read(fasta_path, "fasta")
    sequence = str(record.seq)
    sequence = sequence.upper()
    ori_sequence = sequence[consensus_start:consensus_end]

    return ori_sequence, consensus_start, consensus_end
