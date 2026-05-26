def reverse_complement(dna_seq):
    complement = str.maketrans("ATGC", "TACG")
    rev_dna = dna_seq.translate(complement)[::-1]
    return rev_dna


print(reverse_complement("ATGCCCTAGGGGATGTTTTAA"))


def reading_frames(dna_seq):
    dna_seq = dna_seq.upper()
    frames = []
    # get reverse complement STRING
    rev_seq = reverse_complement(dna_seq)
    # forward frames
    # “Start at position 0, take 3 letters at a time, and keep jumping by 3 until the end
    frames.append([dna_seq[i : i + 3] for i in range(0, len(dna_seq), 3)])
    frames.append([dna_seq[i : i + 3] for i in range(1, len(dna_seq), 3)])
    frames.append([dna_seq[i : i + 3] for i in range(2, len(dna_seq), 3)])
    # reverse frames
    frames.append([rev_seq[i : i + 3] for i in range(0, len(rev_seq), 3)])
    frames.append([rev_seq[i : i + 3] for i in range(1, len(rev_seq), 3)])
    frames.append([rev_seq[i : i + 3] for i in range(2, len(rev_seq), 3)])
    return frames


print(reading_frames("ATGCCCTAGGGGATGTTTTAA"))


# Loop the 6 open reading frames and extract ORFs.
def find_orfs(frames):
    stop_codons = {"TAA", "TAG", "TGA"}
    all_orfs = []

    for frame_index, frame in enumerate(frames):
        i = 0
        while i < len(frame):
            if frame[i] == "ATG":
                orf = ["ATG"]
                j = i + 1
                while j < len(frame):
                    codon = frame[j]
                    orf.append(codon)
                    if codon in stop_codons:
                        all_orfs.append((frame_index + 1, "".join(orf)))
                        break
                    j += 1
                i = j  # move forward after stop
            else:
                i += 1
    return all_orfs


# assuming reverse_complement is defined
frames = reading_frames("ATGCCCTAGGGGATGTTTTAA")
orfs = find_orfs(frames)
