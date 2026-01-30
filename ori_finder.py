from Bio import SeqIO
import statistics


def compute_gc_skew(fragment):
    g = fragment.count("G")
    c = fragment.count("C")
    if g + c == 0:
        return 0
    return (g - c) / (g + c)


def scan_ori_region(fasta_path, window_size=200, step_size=50):
    record = SeqIO.read(fasta_path, "fasta")
    sequence = str(record.seq).upper()

    skew_profile = []
    coord_list = []

    for idx in range(0, len(sequence) - window_size, step_size):
        chunk = sequence[idx:idx + window_size]
        skew_profile.append(compute_gc_skew(chunk))
        coord_list.append(idx)

    change = []
    for i in range(len(skew_profile) - 1):
        change.append(abs(skew_profile[i + 1] - skew_profile[i]))

    peak_index = change.index(max(change))
    start = coord_list[peak_index]
    end = start + window_size

    return start, end


def find_ori_multi_scale(fasta_path, windows=(150, 200, 250, 300), step=25):
    start_points = []
    end_points = []

    for w in windows:
        s, e = scan_ori_region(fasta_path, window_size=w, step_size=step)
        start_points.append(s)
        end_points.append(e)

    consensus_start = int(statistics.median(start_points))
    consensus_end = int(statistics.median(end_points))

    record = SeqIO.read(fasta_path, "fasta")
    sequence = str(record.seq).upper()
    ori_sequence = sequence[consensus_start:consensus_end]

    return ori_sequence, consensus_start, consensus_end
