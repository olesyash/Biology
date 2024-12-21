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



def find_nonsense(position):
    res = 0
    stop_codon = [key for key, val in gencode.items() if val == '_']
    for key, val in gencode.items():
        if val == '_':
            continue
        letters = ['A', 'T', 'C', 'G']
        for l in letters:
            new_key = key[:position-1] + l + key[position:]
            if new_key in stop_codon:
                res += 1
    return res

def find_synonymous(position):
    res = 0
    for key, val in gencode.items():
        letters = ['A', 'T', 'C', 'G']
        for l in letters:
            new_key = key[:position-1] + l + key[position:]
            if new_key == key:
                continue
            if gencode[new_key] == val:
                res += 1
    return res


def count_mutation_by_type(position,type):
    if type == 'nonsense':
        return find_nonsense(position=position)
    if type == 'synonymous':
        return find_synonymous(position)


print(count_mutation_by_type(1, 'synonymous'))