import operator

from Bio import SeqIO
from vaxrank.manufacturability import ManufacturabilityScores


def get_all_kmers(seq, k):
    out = []
    for i in range(len(seq) - k + 1):
        out.append(str(seq[i:i+k]))
    return out


def get_all_kmers_with_overlap(seq, k, overlap):
    out = []
    for i in range(0, len(seq) - k + 1, overlap):
        out.append(str(seq[i:i+k]))
    # this doesn't divide evenly, get the last kmer
    out.append(str(seq[-k:]))
    return out


def print_subsequences(rbd_region):
    print('\n50mers, 25aa overlap')
    for x in get_all_kmers_with_overlap(rbd_region, 50, 25):
        print(x)

    # get all 40mers with 20aa overlap from RBD
    print('\n40mers, 20aa overlap')
    for x in get_all_kmers_with_overlap(rbd_region, 40, 20):
        print(x)

    print('\n30mers, 15aa overlap')
    for x in get_all_kmers_with_overlap(rbd_region, 30, 15):
        print(x)

    print('\n20mers, 10aa overlap')
    for x in get_all_kmers_with_overlap(rbd_region, 20, 10):
        print(x)


def save_21mers(rbd_region):
    all_21mers = get_all_kmers(rbd_region, 21)
    print('%d 21mers' % len(all_21mers))

    # write file, one 21mer per line
    with open('rbd_21mers.txt', 'w') as f:
        for kmer in all_21mers:
            f.write(kmer + '\n')


def main(args_list=None):

    fasta_path = 'spike-glycoprotein-6vsb.fasta'
    seq_records = list(SeqIO.parse(fasta_path, 'fasta'))

    full_seq = seq_records[0].seq
    rbd_seq = seq_records[1].seq

    print('Full seq')
    print(full_seq)
    
    print('RBD seq')
    print(rbd_seq)
    rbd_start = full_seq.find(rbd_seq)

    # get all 21mers encompassing the region [RBD start - 50, RBD end + 50], save to file
    flank_left = rbd_start - 50
    flank_right = rbd_start + len(rbd_seq) + 50
    rbd_region = full_seq[flank_left:flank_right]
    print('RBD +/- flank')
    print(rbd_region)
    print('\n')


    
    # save_21mers(rbd_region)
    # print_subsequences(rbd_region)

    # rank by manufacturability scores
    # seq_to_score = {}
    # for seq in all_21mers:
    #     seq_to_score[seq] = ManufacturabilityScores.from_amino_acids(seq)
    # sorted_seqs = sorted(seq_to_score.items(), key=operator.itemgetter(1))
    


if __name__ == "__main__":
    main()
