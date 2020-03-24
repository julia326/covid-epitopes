import operator

from Bio import SeqIO
from vaxrank.manufacturability import ManufacturabilityScores


def get_all_kmers(seq, k):
    out = []
    for i in range(len(seq) - k + 1):
        out.append(str(seq[i:i+k]))
    return out


def main(args_list=None):

    fasta_path = 'spike-glycoprotein-6vsb.fasta'
    seq_records = list(SeqIO.parse(fasta_path, 'fasta'))

    full_seq = seq_records[0].seq
    rbd_seq = seq_records[1].seq
    
    print('RBD seq: %s' % rbd_seq)
    rbd_start = full_seq.find(rbd_seq)

    # get all 21mers encompassing the region [RBD start - 50, RBD end + 50], save to file
    flank_left = rbd_start - 50
    flank_right = rbd_start + len(rbd_seq) + 50
    rbd_region = full_seq[flank_left:flank_right]

    all_21mers = get_all_kmers(rbd_region, 21)
    print('%d 21mers' % len(all_21mers))

    # write file, one 21mer per line
    with open('rbd_21mers.txt', 'w') as f:
        for kmer in all_21mers:
            f.write(kmer + '\n')

    # rank by manufacturability scores
    # seq_to_score = {}
    # for seq in all_21mers:
    #     seq_to_score[seq] = ManufacturabilityScores.from_amino_acids(seq)
    # sorted_seqs = sorted(seq_to_score.items(), key=operator.itemgetter(1))
    


if __name__ == "__main__":
    main()
