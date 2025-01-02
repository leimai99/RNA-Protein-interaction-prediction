import sys,re
from collections import Counter
import pandas as pd
import itertools
import numpy as np
ALPHABET='ACGU'

def readRNAFasta(file):
	with open(file) as f:
		records = f.read()

	if re.search('>', records) == None:
		print('The input RNA sequence must be fasta format.')
		sys.exit(1)
	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ACGU-]', '-', ''.join(array[1:]).upper())
		myFasta.append([name, sequence])
	return myFasta

def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer
def kmer(input_data, k=2, normalize=True):
    
    fastas=readRNAFasta(input_data)
    vector = []
    header = ['#']
    if k < 1:
        print('error, the k must be positive integer.')
        return 0
    for kmer in itertools.product(ALPHABET, repeat=k):
        header.append(''.join(kmer))
    vector.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        kmers = kmerArray(sequence, k)
        count = Counter()
        count.update(kmers)
        if normalize == True:
           for key in count:
               count[key] = count[key] / len(kmers)
        code = [name]
        for j in range(1, len(header)):
            if header[j] in count:
               code.append(count[header[j]])
            else:
                code.append(0)
        vector.append(code)
    return vector
def read_fasta_file(fasta_file):
    seq_dict = {}
    fp = open(fasta_file, 'r')
    name = ''
    # pdb.set_trace()
    for line in fp:
        # let's discard the newline at the end (if any)
        line = line.rstrip().strip('*')
        # distinguish header from sequence
        if line[0] == '>':  # or line.startswith('>')
            # it is the header
            name = line[1:] # discarding the initial >
            seq_dict[name] = ''
        else:
            # it is sequence
            seq_dict[name] = seq_dict[name] + line
    fp.close()
    for key,value in seq_dict.items():
        if len(value)>10000:
            print(key)
            print(len(value))
    print(len(seq_dict))
    return seq_dict

def get_4_nucleotide_composition(tris, seq, pythoncount=True):
 
    seq_len = len(seq)
    tri_feature = []

    if pythoncount:
        for val in tris:
            num = seq.count(val)
            tri_feature.append(float(num) / seq_len)
    else:
        k = len(tris[0])
        tmp_fea = [0] * len(tris)
        for x in range(len(seq) + 1 - k):
            kmer = seq[x:x + k]
            if kmer in tris:
                ind = tris.index(kmer)
                tmp_fea[ind] = tmp_fea[ind] + 1
        tri_feature = [float(val) / (len(seq) + 1 - k) for val in tmp_fea] 
        # pdb.set_trace()
    return tri_feature

def TransDict_from_list(groups):

    tar_list = ['0', '1', '2', '3']
    result = {}
    index = 0
    # groups = ['AGV', 'ILFP', 'YMTS', 'HNQW', 'RK', 'DE', 'C']
    for group in groups:
        g_members = sorted(group)
        print(g_members)
        for c in g_members:
            result[c] = str(tar_list[index])
        index = index + 1

    return result

def translate_sequence (seq, TranslationDict):

    from_list = []
    to_list = []
    for k,v in TranslationDict.items():
        from_list.append(k)
        to_list.append(v)
    TRANS_seq = seq.translate(str.maketrans(str(from_list), str(to_list)))

    return TRANS_seq

def find_all_path(path,cnt,base,all_kmers,k):
    if(cnt>k):
        all_kmers.append(path)
    else:
        for i in range(base):
            path+=str(i)
            find_all_path(path,cnt+1,base,all_kmers,k)
            path = path[:-1]

def get_k_RNA_trids(k):
  
    chars = ['0', '1', '2', '3']
    base = len(chars)
    all_kmers = []
    cnt = 1
    for j in range(base):
        path = ""
        path += str(j)
        find_all_path(path, cnt + 1, base, all_kmers,k)

    return all_kmers


def generated_RNA_kmer(fasta_file, savepath, k = 4):
    RNA4mer = []
    RNA_seq_dict = read_fasta_file(fasta_file)
    groups = ['A', 'C', 'G', 'U']
    group_dict = TransDict_from_list(groups)
    RNA_tris = get_k_RNA_trids(4)

    for RNA, RNA_seq in RNA_seq_dict.items():
        RNA_seq1 = translate_sequence(RNA_seq, group_dict)
        RNA_tri_fea = get_4_nucleotide_composition(RNA_tris, RNA_seq1, pythoncount=False)
        RNA4mer.append(RNA_tri_fea)
    RNA4mer = np.array(RNA4mer)

    print(RNA4mer.shape)
    RNA4mer = pd.DataFrame(RNA4mer)
    #RNA4mer.to_csv(savepath, index=False)


if __name__ == '__main__':

    generated_RNA_kmer(
        fasta_file = 'D:/CSU/RPI-MD-/RPI-MD/data/generated_data/NPInter_10412/RNA_extracted_seq.fasta',
        savepath = 'D:/CSU/RPI-MD-/RPI-MD/data/generated_data/NPInter_10412/RNA4merfeat.csv',k=4)
    generated_RNA_kmer(
        fasta_file='D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI369/Protein_extracted_seq.fasta',
        savepath='D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI369/RNA4merfeat.csv', k=4)
    generated_RNA_kmer(fasta_file = 'D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI2241/RNA_extracted_seq.fasta',
        savepath = 'D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI2241/RNA4merfeat.csv',k=4)
    generated_RNA_kmer(fasta_file = 'D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI7317/RNA_extracted_seq.fasta',
        savepath = 'D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI7317/RNA4merfeat.csv',k=4)
    generated_RNA_kmer(
        fasta_file='D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI1807/RNA_extracted_seq.fasta',
        savepath='D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI1807/RNA4merfeat.csv', k=4)
kmer_vector=kmer('D:/CSU/RPI-MD-/RPI-MD/data/raw_data/RPI2241_rna_seq.fa',k=4)
csv_data=pd.DataFrame(data=kmer_vector)
csv_data.to_csv('D:/CSU/RPI-MD-/RPI-MD/data/raw_data/RPI2241.csv')