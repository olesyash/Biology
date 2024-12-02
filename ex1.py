# Olesya Sharify, 319346565
# Adi Aharoni, 211749361

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}


# q1: Check if the seq is a valid DNA
def is_valid_dna(seq):
    legal_letters = {'t', 'g', 'c', 'a', 'T', 'G', 'C', 'A'}
    for letter in seq:
        if letter not in legal_letters:
            return False
    return True


# q2: Count the appearance of the letters G and C in the seq
# returns the percentage of their appearance
def get_gc_content(seq):
    if not is_valid_dna(seq):
        return -1
    count = 0
    count_vals = {'C', 'c', 'G', 'g'}
    for letter in seq:
        if letter in count_vals:
            count += 1
    return round(100*(count/len(seq)), 2)


# q3: Returns the reverse complement of a given seq
# 1. reverse the seq
# 2. turning the letters: A -> T, T -> A, C -> G, G -> C
def reverse_complement(seq):
    if not is_valid_dna(seq):
        return ''
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_seq = seq[::-1]
    res = ''
    for letter in rev_seq:
        res += complement[letter.upper()]
    return res


# q4: turning a seq into RNA
# * If seq = 1: all we ave to do is change T to U (we're on the positive strand).
# * If seq = -1: reverse and thank change T to U (we're on the negative strand).
# * We assume strand value is legal; 1/-1.
def get_transcription(seq, strand):
    if not is_valid_dna(seq):
        return ''
    if strand == 1:
        return seq.upper().replace('T', 'U')
    else:
        return reverse_complement(seq).upper().replace('T', 'U')


# q5: count the appearance of each amino in the given gencode
def stats_amino_acids():
    global gencode
    res = {}
    gen_vals = list(gencode.values())
    for unique_letter in set(gen_vals):
        res[unique_letter] = gen_vals.count(unique_letter)
    if '_' in res:
        res['END'] = res.pop('_')
    return res


# q6: Translating a given seq into a sequence of amino acids
# * starting from letter number {reading_frame = 1/2/3} in the original seq.
# * is reading_frame = None: returning an array with all the different possibilities of the translation.
# * We assume reading_frame value will bw valid: 1/2/3/None.
def translate_seq(seq, reading_frame=None):
    if not is_valid_dna(seq):
        return []
    global gencode
    res = []
    if not reading_frame:
        # Translate all three possible reading frames
        for i in range(3):
            protein = ''
            for j in range(i, len(seq) - 2, 3):
                cur = seq[j:j+3]
                amino = gencode[cur]
                if amino.upper() == '_':
                    break
                protein += amino
            res.append(protein)
        return res
    start = reading_frame - 1
    protein = ''
    for j in range(start, len(seq) - 2, 3):
        cur = seq[j:j+3]
        amino = gencode[cur]
        if amino.upper() == '_':
            break
        protein += amino
    res.append(protein)
    return res


# q7: count the amino acids appearance in the seq
# * only the amino acids that appear in the seq will be shown
# * the reading starts from the first letter in the seq.
def count_codons(seq):
    if not is_valid_dna(seq):
        return {}
    res = {}
    for i in range(0, len(seq) - 2, 3):
        cur = seq[i:i+3].upper()
        if cur in res:
            res[cur] += 1
        else:
            res[cur] = 1
    return res


def main():
    # d1:
    print('Demonstration the function on a valid and invalid seq:')
    seq1 = 'aaTgTGAtcCc'
    print(f'seq1 = {seq1}')
    print(f'seq1 is valid dna? {is_valid_dna(seq1)}')
    seq2 = 'aaTgpGAtcCc'
    print(f'seq2 = {seq1}')
    print(f'seq2 is valid dna? {is_valid_dna(seq2)}')

    # d2:
    print('''\ndemonstration of get_gc_content function:
     1. On the valid seq1 = aaTgTGAtcCc.
     2. On the invalid seq2 = aaTgpGAtcCc. ''')
    print(f'seq1 gc content = {get_gc_content(seq1)}%')
    seq2_gc_content = get_gc_content(seq2)
    if seq2_gc_content == -1:
        print("seq2 is not valid DNA sequence")
    else:
        print(f'seq2 gc content = {seq2_gc_content}%')

    # d3:
    print('\nDemonstration of reverse_complement function:')
    seq1 = 'ACgTtG'
    print(f'the reverse complement of the seq: ACgTtG is: {reverse_complement(seq1)}')

    # d4: Demonstrating the RNA of a sequence while:
    # 1. We are on the positive strand
    # 2. We are on the negative strand
    print('\nDemonstrating the RNA of the sequence: ACgTtGa')
    seq1 = 'ACgTtGa'
    print(f'positive strand example: {get_transcription(seq1, 1)}')
    print(f'negative strand example: {get_transcription(seq1, -1)}')

    # d5:
    print('\nDemonstrating the result of stats_amino_acids on given gencode:')
    print(stats_amino_acids())

    # d6:
    print('''\nDemonstration of translate_seq function on the seq: ATACCCTAGTGT: 
    1. Start reading from the first letter of the seq.
    2. Start reading from the second letter of the seq.
    3. Start reading from the third letter of the seq.
    4. Array of all the different possibilities described above.''')
    print(f'The given gencode is: {gencode}')
    print(translate_seq('ATACCCTAGTGT', 1))
    print(translate_seq('ATACCCTAGTGT', 2))
    print(translate_seq('ATACCCTAGTGT', 3))
    print(translate_seq('ATACCCTAGTGT'))

    # d7:
    print('\nDemonstration of count_codons function on the seq: ACGTGTTGTTTcAAT')
    print(count_codons('ACGTGTTGTTTcAAT'))


if __name__ == "__main__":
    main()