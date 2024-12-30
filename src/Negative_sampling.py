import numpy as np
import pandas as pd
import math
import random

np.random.seed(1)
random.seed(1)

def calculate_protein_sw_similarity(pr1, pr2, swscore_matrix):

    score = swscore_matrix[pr1, pr2]/math.sqrt(swscore_matrix[pr1, pr1]*swscore_matrix[pr2, pr2])
    return score

def calculate_socre_of_pri_and_RNAj(pr_i, RNA_j,positive_samples,swscore_matrix):

    score = 0
    related_pair = [pair for pair in positive_samples if pair[0] == RNA_j]
    for pair in related_pair:
        if(pair[1]!= pr_i):
            score += calculate_protein_sw_similarity(pr_i,pair[1],swscore_matrix)
    return score


def get_positive_samples_of_NPInter(NPI_filepath):
    NPInter = pd.read_table(NPI_filepath)
    protein = NPInter['UNIPROT-ID'].unique().tolist()  
    RNA = NPInter['NONCODE-ID'].unique().tolist()  
    positive_index = []  
    for index, row in NPInter.iterrows():
        i = RNA.index(row['NONCODE-ID'])
        j = protein.index(row['UNIPROT-ID'])
        positive_index.append([i, j])
    return positive_index, protein, RNA


def get_Positives_and_Negatives (positive_samples, pr_list, RNA_list,swscore_matrix, savepath):
    Positives = []
    Negatives = []

    # random pair
    for RNA_index in range((len(RNA_list))):
        for pr_index in range(len(pr_list)):
            sample = [RNA_index, pr_index]
            if [RNA_index, pr_index] in positive_samples:
                Ms = 1
                sample.append(Ms)
                Positives.append(sample)
            else:
                Ms = calculate_socre_of_pri_and_RNAj(pr_index, RNA_index, positive_samples,swscore_matrix)
                sample.append(Ms)
                Negatives.append(sample)
    Negatives = sorted(Negatives, key=lambda x: x[2])

    Positives = pd.DataFrame(Positives, columns=['RNA', 'protein', 'label'])
    Negatives = pd.DataFrame(Negatives, columns=['RNA', 'protein', 'label'])
    Positives.to_csv(savepath+'Positives.csv', index=False)
    Negatives.to_csv(savepath+'Negatives.csv', index=False)

    return Positives, Negatives

def get_edgelist(Positives, Negatives, method, savepath, nc_num):

    if method == 'sort':
        Negatives = Negatives[:(len(Positives))]  
    elif method == 'random':
   #     Negatives = Negatives.sample(n=len(Positives), random_state=1) 
   # elif method == 'sort_random':

        Negatives = Negatives[:len(Positives) * 2]
        Negatives = Negatives.sample(n=len(Positives), random_state=1)
    elif method =='raw':

        pass
    Negatives.loc[:, 'label'] = -1

    edgelist = [Positives, Negatives]
    edgelist = pd.concat(edgelist, axis=0)
    edgelist = edgelist.take(np.random.permutation(len(edgelist)))
    edgelist = edgelist.reset_index(drop=True)

    if method == 'sort':
        edgelist.to_csv(savepath + 'edgelist_sort.csv', header=None)
        edgelist['protein'] = edgelist['protein'] + nc_num
        edgelist = edgelist[['RNA', 'protein']]
        np.savetxt(savepath + 'graph.edgelist_sort.txt', edgelist, fmt='%s',
               delimiter=' ')
    elif method == 'random':
        edgelist.to_csv(savepath + 'edgelist_random.csv', header=None)
        edgelist['protein'] = edgelist['protein'] + nc_num
        edgelist = edgelist[['RNA', 'protein']]
        np.savetxt(savepath + 'graph.edgelist_random.txt', edgelist, fmt='%s',
                   delimiter=' ')
    elif method == 'sort_random':
        edgelist.to_csv(savepath + 'edgelist_sort_random.csv', header=None)
        edgelist['protein'] = edgelist['protein'] + nc_num
        edgelist = edgelist[['RNA', 'protein']]
        np.savetxt(savepath + 'graph.edgelist_sort_random.txt', edgelist, fmt='%s',
                   delimiter=' ')
    elif method == 'raw':
        edgelist.to_csv(savepath + 'edgelist_raw.csv', header=None)
        edgelist['protein'] = edgelist['protein'] + nc_num
        edgelist = edgelist[['RNA', 'protein']]
        np.savetxt(savepath + 'graph.edgelist_raw.txt', edgelist, fmt='%s',
                   delimiter=' ')

    return Positives, Negatives


def get_NPInter(filepath, savepath):

    swpath = savepath + 'protein sw_smilarity matrix.csv'
    swscore_matrix = pd.read_csv(swpath, header=None).values
    positive_samples, protein, RNA = get_positive_samples_of_NPInter(filepath)
    #Positives,Negatives = get_Positives_and_Negatives(positive_samples, protein, RNA, swscore_matrix, savepath)
    Positives = pd.read_csv(savepath + 'Positives.csv')
    Negatives = pd.read_csv(savepath + 'Negatives.csv')  

    Positives, Negatives_sort = get_edgelist(Positives, Negatives,'sort',savepath, len(RNA)) 
    _, Negatives_random = get_edgelist(Positives, Negatives,'random', savepath, len(RNA))
    _, Negatives_sort_random = get_edgelist(Positives, Negatives, 'sort_random', savepath, len(RNA))

    NPI_pos = np.zeros((len(RNA), len(protein)))
    NPI_pos[Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    NPI_neg_sort, NPI_neg_random ,NPI_neg_sort_random  =  np.zeros((len(RNA), len(protein))),  np.zeros((len(RNA), len(protein))), np.zeros((len(RNA), len(protein)))
    NPI_neg_sort[Negatives_sort.values[:, 0], Negatives_sort.values[:, 1]] = 1
    NPI_neg_sort = pd.DataFrame(NPI_neg_sort)
    NPI_neg_sort.to_csv(savepath+'NPI_neg_sort.csv', index=False, header=None)

    NPI_neg_random[Negatives_random.values[:, 0], Negatives_random.values[:, 1]] = 1
    NPI_neg_random = pd.DataFrame(NPI_neg_random)
    NPI_neg_random.to_csv(savepath + 'NPI_neg_random.csv', index=False, header=None)

    NPI_neg_sort_random[Negatives_sort_random.values[:, 0], Negatives_sort_random.values[:, 1]] = 1
    NPI_neg_sort_random = pd.DataFrame(NPI_neg_sort_random)
    NPI_neg_sort_random.to_csv(savepath + 'NPI_neg_sort_random.csv', index=False, header=None)


def get_RPI13254(filepath, savepath):

    RPI13254_positive = pd.read_table(filepath + 'RPI13254_positive.txt')
    RPI13254_negative = pd.read_table(filepath + 'RPI13254_negative.txt')

    RPI13254_pos = RPI13254_positive['gene'].tolist()
    protein_pos = [item[:7] for item in RPI13254_pos] 
    RNA_pos = [item[8:] for item in RPI13254_pos]  

    RPI13254_neg = RPI13254_negative['gene'].tolist()
    protein_neg = [item[:7] for item in RPI13254_neg] 
    RNA_neg = [item[8:] for item in RPI13254_neg]  

    protein = list(set(protein_pos).union(set(protein_neg))) 
    RNA = list(set(RNA_pos).union(set(RNA_neg)))  

    merge1 = {"RNA": RNA_pos,"protein": protein_pos}
    RPI13254_pos = pd.DataFrame(merge1)  
    merge2 = {"RNA": RNA_neg, "protein": protein_neg}
    RPI13254_neg = pd.DataFrame(merge2)

    discard_list = ['YBL039W-A', 'YBL101W-A', 'YFL057C', 'YIR044C', 'YAR062W', 'YNL097C-A']
    index_list = RPI13254_pos[RPI13254_pos['RNA'].isin(discard_list)].index
    RPI13254_pos = RPI13254_pos.drop(index=index_list)
    RNA = RPI13254_pos['RNA'].unique().tolist()
    protein = RPI13254_pos['protein'].unique().tolist()
    print("RPI13254:")
    print("protein:" + str(len(protein)) + " RNA:" + str(len(RNA)) + " positives:" + str(
        len(RPI13254_pos)) + " negatives:" + str(len(RPI13254_neg)))

    print(RNA)
    print(protein)

    positive_index = []
    negative_index = []

    for index, row in RPI13254_pos.iterrows():
        i = RNA.index(row['RNA'])
        j = protein.index(row['protein'])
        positive_index.append([i, j])

    for index, row in RPI13254_neg.iterrows():
        i = RNA.index(row['RNA'])
        j = protein.index(row['protein'])
        negative_index.append([i, j, -1])

    Negatives_raw = pd.DataFrame(negative_index,columns = ['RNA', 'protein', 'label'])
    Negatives_raw.to_csv(savepath+'Negatives_raw.csv', index=False)

    swpath = savepath + 'protein sw_smilarity matrix.csv'
    swscore_matrix = pd.read_csv(swpath, header=None).values
    Positives, Negatives = get_Positives_and_Negatives(positive_index, protein, RNA, swscore_matrix, savepath)
    Positives = pd.read_csv(savepath + 'Positives.csv')
    Negatives = pd.read_csv(savepath + 'Negatives.csv')

    Positives, Negatives_sort = get_edgelist(Positives, Negatives, 'sort', savepath,
                                             len(RNA)) 
    _, Negatives_raw = get_edgelist(Positives, Negatives_raw, 'raw', savepath, len(RNA))
    _, Negatives_random = get_edgelist(Positives, Negatives, 'random', savepath, len(RNA))
    _, Negatives_sort_random = get_edgelist(Positives, Negatives, 'sort_random', savepath, len(RNA))

    NPI_pos = np.zeros((len(RNA), len(protein)))
    NPI_pos[Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    NPI_neg_sort, NPI_neg_raw, NPI_neg_random, NPI_neg_sort_random = np.zeros((len(RNA), len(protein))), np.zeros(
        (len(RNA), len(protein))), np.zeros((len(RNA), len(protein))), np.zeros((len(RNA), len(protein)))

    NPI_neg_sort[Negatives_sort.values[:, 0], Negatives_sort.values[:, 1]] = 1
    NPI_neg_sort = pd.DataFrame(NPI_neg_sort)
    NPI_neg_sort.to_csv(savepath + 'NPI_neg_sort.csv', index=False, header=None)

    NPI_neg_raw[Negatives_raw.values[:, 0], Negatives_raw.values[:, 1]] = 1
    NPI_neg_raw = pd.DataFrame(NPI_neg_raw)
    NPI_neg_raw.to_csv(savepath + 'NPI_neg_raw.csv', index=False, header=None)

    NPI_neg_random[Negatives_random.values[:, 0], Negatives_random.values[:, 1]] = 1
    NPI_neg_random = pd.DataFrame(NPI_neg_random)
    NPI_neg_random.to_csv(savepath + 'NPI_neg_random.csv', index=False, header=None)

    NPI_neg_sort_random[Negatives_sort_random.values[:, 0], Negatives_sort_random.values[:, 1]] = 1
    NPI_neg_sort_random = pd.DataFrame(NPI_neg_sort_random)
    NPI_neg_sort_random.to_csv(savepath + 'NPI_neg_sort_random.csv', index=False, header=None)


def get_RPI7317(filepath, savepath):
    RPI7317 = pd.read_csv(filepath)
    protein = RPI7317['Protein names'].unique().tolist()
    RNA = RPI7317['RNA names'].unique().tolist()
    positive_index = []
    print(protein)
    print(RNA)
    for index, row in RPI7317.iterrows():
        i = RNA.index(row['RNA names'])
        j = protein.index(row['Protein names'])
        positive_index.append([i, j])

    print("RPI7317")
    print("protein:" + str(len(protein)) + " RNA:" + str(len(RNA)) + " positives:" + str(len(positive_index)))
    swpath = savepath + 'protein sw_smilarity matrix.csv'
    swscore_matrix = pd.read_csv(swpath, header=None).values

    Positives = pd.read_csv(savepath + 'Positives.csv')
    Negatives = pd.read_csv(savepath + 'Negatives.csv')

    Positives, Negatives_sort = get_edgelist(Positives, Negatives, 'sort', savepath,
                                             len(RNA))  
    _, Negatives_random = get_edgelist(Positives, Negatives, 'random', savepath, len(RNA))
    _, Negatives_sort_random = get_edgelist(Positives, Negatives, 'sort_random', savepath, len(RNA))

    NPI_pos = np.zeros((len(RNA), len(protein)))
    NPI_pos[Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    NPI_neg_sort, NPI_neg_random, NPI_neg_sort_random = np.zeros((len(RNA), len(protein))), np.zeros(
        (len(RNA), len(protein))), np.zeros((len(RNA), len(protein)))
    NPI_neg_sort[Negatives_sort.values[:, 0], Negatives_sort.values[:, 1]] = 1
    NPI_neg_sort = pd.DataFrame(NPI_neg_sort)
    NPI_neg_sort.to_csv(savepath + 'NPI_neg_sort.csv', index=False, header=None)

    NPI_neg_random[Negatives_random.values[:, 0], Negatives_random.values[:, 1]] = 1
    NPI_neg_random = pd.DataFrame(NPI_neg_random)
    NPI_neg_random.to_csv(savepath + 'NPI_neg_random.csv', index=False, header=None)

    NPI_neg_sort_random[Negatives_sort_random.values[:, 0], Negatives_sort_random.values[:, 1]] = 1
    NPI_neg_sort_random = pd.DataFrame(NPI_neg_sort_random)
    NPI_neg_sort_random.to_csv(savepath + 'NPI_neg_sort_random.csv', index=False, header=None)


def get_RPI1446(filepath, savepath):
    RPI1446 = pd.read_table(filepath, header=None, names=['protein', 'RNA','label'])
    protein = RPI1446['protein'].unique().tolist()
    RNA = RPI1446['RNA'].unique().tolist()
    positive_index = []
    negative_index = []

    for index,row in RPI1446.iterrows():
        i = RNA.index(row['RNA'])
        j = protein.index(row['protein'])
        if(row['label']==1):
            positive_index.append([i,j])
        else:
            negative_index.append([i,j,-1])
    print("RPI1446:")
    print("positive:"+str(len(positive_index))+" negative:"+str(len(negative_index))+" RNA:"+str(len(RNA))+" protein:"+str(len(protein)))
    Negatives_raw = pd.DataFrame(negative_index, columns=['RNA', 'protein', 'label'])
    Negatives_raw.to_csv(savepath + 'Negatives_raw.csv', index=False)
    Positives = pd.DataFrame(positive_index, columns =['RNA', 'protein'])
    Positives['label'] = 1
    get_edgelist(Positives, Negatives_raw, 'raw', savepath, len(RNA))

    NPI_pos = np.zeros((len(RNA), len(protein)))
    NPI_pos[Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    NPI_neg_raw = np.zeros((len(RNA), len(protein)))
    NPI_neg_raw[Negatives_raw.values[:, 0], Negatives_raw.values[:, 1]] = 1
    NPI_neg_raw = pd.DataFrame(NPI_neg_raw)
    NPI_neg_raw.to_csv(savepath + 'NPI_neg_raw.csv', index=False, header=None)


def get_RPI1807(filepath, savepath):
    RPI1807_positive = pd.read_table(filepath + 'RPI1807_PositivePairs.csv')
    RPI1807_negative = pd.read_table(filepath + 'RPI1807_NegativePairs.csv')
    protein = list(set(RPI1807_positive['Protein ID'].tolist()).union(set(RPI1807_negative['Protein ID'].tolist())))
    RNA = list(set(RPI1807_positive['RNA ID'].tolist()).union(set(RPI1807_negative['RNA ID'].tolist())))
    protein_pos = RPI1807_positive['Protein ID'].unique().tolist()
    RNA_pos = RPI1807_positive['RNA ID'].unique().tolist()
    positive_index = []  
    negative_index = []

    for index, row in RPI1807_positive.iterrows():
        i = RNA.index(row['RNA ID'])
        j = protein.index(row['Protein ID'])
        positive_index.append([i, j])

    for index, row in RPI1807_negative.iterrows():
        i = RNA.index(row['RNA ID'])
        j = protein.index(row['Protein ID'])
        negative_index.append([i, j, -1])

    Positives = pd.DataFrame(positive_index, columns = ['RNA','protein'])
    Positives.loc[:,'label'] = 1
    Negatives_raw = pd.DataFrame(negative_index,columns = ['RNA','protein','label'])
    Positives.to_csv(savepath +'Positives.csv',index=False)
    Negatives_raw.to_csv(savepath + 'Negatives_raw.csv',index=False)
    _, Negatives_raw = get_edgelist(Positives, Negatives_raw, 'raw', savepath, len(RNA))
    print("RPI1807")
    print("protein:" + str(len(protein)) + " RNA:" + str(len(RNA)) + " positives:" + str(len(positive_index)))

    positive_index = []
    for index, row in RPI1807_positive.iterrows():
        i = RNA_pos.index(row['RNA ID'])
        j = protein_pos.index(row['Protein ID'])
        positive_index.append([i, j])

    swpath = savepath + 'protein sw_smilarity matrix.csv'
    swscore_matrix = pd.read_csv(swpath, header=None).values
    Positives, Negatives = get_Positives_and_Negatives(positive_index, protein, RNA, swscore_matrix, savepath)

    Positives, Negatives_sort = get_edgelist(Positives, Negatives, 'sort', savepath,
                                             len(RNA)) 
    _, Negatives_random = get_edgelist(Positives, Negatives, 'random', savepath, len(RNA))
    _, Negatives_sort_random = get_edgelist(Positives, Negatives, 'sort_random', savepath, len(RNA))

    NPI_pos, NPI_neg_random, NPI_neg_sort, NPI_neg_sort_random, NPI_neg_raw = np.zeros((len(RNA), len(protein))), np.zeros((len(RNA), len(protein))),\
                                                                              np.zeros((len(RNA), len(protein))), np.zeros((len(RNA), len(protein))), np.zeros((len(RNA), len(protein)))
    NPI_pos [Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    NPI_neg_sort[Negatives_sort.values[:, 0], Negatives_sort.values[:, 1]] = 1
    NPI_neg_sort = pd.DataFrame(NPI_neg_sort)
    NPI_neg_sort.to_csv(savepath + 'NPI_neg_sort.csv', index=False, header=None)

    NPI_neg_raw[Negatives_raw.values[:, 0], Negatives_raw.values[:, 1]] = 1
    NPI_neg_raw = pd.DataFrame(NPI_neg_raw)
    NPI_neg_raw.to_csv(savepath + 'NPI_neg_raw.csv', index=False, header=None)

    NPI_neg_random[Negatives_random.values[:, 0], Negatives_random.values[:, 1]] = 1
    NPI_neg_random = pd.DataFrame(NPI_neg_random)
    NPI_neg_random.to_csv(savepath + 'NPI_neg_random.csv', index=False, header=None)

    NPI_neg_sort_random[Negatives_sort_random.values[:, 0], Negatives_sort_random.values[:, 1]] = 1
    NPI_neg_sort_random = pd.DataFrame(NPI_neg_sort_random)
    NPI_neg_sort_random.to_csv(savepath + 'NPI_neg_sort_random.csv', index=False, header=None)


def get_RPI369(filepath, savepath):
    RPI369 = pd.read_table(filepath , header=None,
                           names=['protein', 'RNA', 'label'])
    RPI369_pos = RPI369[RPI369['label']==1]
    protein = RPI369_pos['protein'].unique().tolist()
    RNA = RPI369_pos['RNA'].unique().tolist()
    positive_index = [] 

    for index, row in RPI369_pos.iterrows():
        i = RNA.index(row['RNA'])
        j = protein.index(row['protein'])
        positive_index.append([i, j])
    print("RPI369:")
    print("positive:" + str(len(positive_index)) + " RNA:" + str(len(RNA)) + " protein:" + str(len(protein)))


    Positives = pd.DataFrame(positive_index, columns=['RNA', 'protein'])
    Positives['label'] = 1

    NPI_pos = np.zeros((len(RNA), len(protein)))
    NPI_pos[Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    swpath = savepath + 'protein sw_smilarity matrix.csv'
    swscore_matrix = pd.read_csv(swpath, header=None).values
    Positives,Negatives = get_Positives_and_Negatives(positive_index, protein, RNA, swscore_matrix, savepath)


    Positives, Negatives_sort = get_edgelist(Positives, Negatives, 'sort', savepath,
                                             len(RNA)) 
    _, Negatives_random = get_edgelist(Positives, Negatives, 'random', savepath, len(RNA))
    _, Negatives_sort_random = get_edgelist(Positives, Negatives, 'sort_random', savepath, len(RNA))

    NPI_neg_sort, NPI_neg_random, NPI_neg_sort_random = np.zeros((len(RNA), len(protein))), np.zeros(
        (len(RNA), len(protein))), np.zeros((len(RNA), len(protein)))
    NPI_neg_sort[Negatives_sort.values[:, 0], Negatives_sort.values[:, 1]] = 1
    NPI_neg_sort = pd.DataFrame(NPI_neg_sort)
    NPI_neg_sort.to_csv(savepath + 'NPI_neg_sort.csv', index=False, header=None)

    NPI_neg_random[Negatives_random.values[:, 0], Negatives_random.values[:, 1]] = 1
    NPI_neg_random = pd.DataFrame(NPI_neg_random)
    NPI_neg_random.to_csv(savepath + 'NPI_neg_random.csv', index=False, header=None)

    NPI_neg_sort_random[Negatives_sort_random.values[:, 0], Negatives_sort_random.values[:, 1]] = 1
    NPI_neg_sort_random = pd.DataFrame(NPI_neg_sort_random)
    NPI_neg_sort_random.to_csv(savepath + 'NPI_neg_sort_random.csv', index=False, header=None)


def get_RPI2241(filepath, savepath):
    RPI2241 = pd.read_table(filepath + 'RPI2241_pairs.txt', header=None, names=['Protein', 'RNA', 'label'])
    protein = RPI2241['Protein'].unique().tolist()
    RNA = RPI2241['RNA'].unique().tolist()
    positive_index = []
    negative_index = []

    for index, row in RPI2241.iterrows():
        i = RNA.index(row['RNA'])
        j = protein.index(row['Protein'])
        if (row['label'] == 1):
            positive_index.append([i, j])
        else:
            negative_index.append([i, j, -1])
    print("RPI2241:")
    print("positive:" + str(len(positive_index)) + " negative:" + str(len(negative_index)) + " RNA:" + str(
        len(RNA)) + " protein:" + str(len(protein)))
    Negatives_raw = pd.DataFrame(negative_index, columns=['RNA', 'protein', 'label'])
    Negatives_raw.to_csv(savepath + 'Negatives_raw.csv', index=False)

    Positives = pd.DataFrame(positive_index, columns=['RNA', 'protein'])
    Positives['label'] = 1
    get_edgelist(Positives, Negatives_raw, 'raw', savepath, len(RNA))

    NPI_pos = np.zeros((len(RNA), len(protein)))
    NPI_pos[Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    NPI_neg_raw = np.zeros((len(RNA), len(protein)))
    NPI_neg_raw[Negatives_raw.values[:, 0], Negatives_raw.values[:, 1]] = 1
    NPI_neg_raw = pd.DataFrame(NPI_neg_raw)
    NPI_neg_raw.to_csv(savepath + 'NPI_neg_raw.csv', index=False, header=None)

    swpath = savepath + 'protein sw_smilarity matrix.csv'
    swscore_matrix = pd.read_csv(swpath, header=None).values
    Positives,Negatives = get_Positives_and_Negatives(positive_index, protein, RNA, swscore_matrix, savepath)


    Positives, Negatives_sort = get_edgelist(Positives, Negatives, 'sort', savepath,
                                             len(RNA)) 
    _, Negatives_random = get_edgelist(Positives, Negatives, 'random', savepath, len(RNA))
    _, Negatives_sort_random = get_edgelist(Positives, Negatives, 'sort_random', savepath, len(RNA))

    NPI_neg_sort, NPI_neg_random, NPI_neg_sort_random = np.zeros((len(RNA), len(protein))), np.zeros(
        (len(RNA), len(protein))), np.zeros((len(RNA), len(protein)))
    NPI_neg_sort[Negatives_sort.values[:, 0], Negatives_sort.values[:, 1]] = 1
    NPI_neg_sort = pd.DataFrame(NPI_neg_sort)
    NPI_neg_sort.to_csv(savepath + 'NPI_neg_sort.csv', index=False, header=None)

    NPI_neg_random[Negatives_random.values[:, 0], Negatives_random.values[:, 1]] = 1
    NPI_neg_random = pd.DataFrame(NPI_neg_random)
    NPI_neg_random.to_csv(savepath + 'NPI_neg_random.csv', index=False, header=None)

    NPI_neg_sort_random[Negatives_sort_random.values[:, 0], Negatives_sort_random.values[:, 1]] = 1
    NPI_neg_sort_random = pd.DataFrame(NPI_neg_sort_random)
    NPI_neg_sort_random.to_csv(savepath + 'NPI_neg_sort_random.csv', index=False, header=None)


def get_NPInter4158(filepath, savepath):
    NPInter4158 = pd.read_table(filepath + 'NPInter4158_interaction.txt', header=None, sep=' ')
    RNA = pd.read_table(filepath + 'NPInter4158_RNA.txt', header=None)
    protein = pd.read_table(filepath + 'NPInter4158_protein.txt', header=None)
    NPInter4158 = NPInter4158.values
    index = NPInter4158.nonzero()
    positive_index = []
    for i, _ in enumerate(index[0]):
        positive_index.append([index[0][i], index[1][i]])

    swpath = savepath + 'protein sw_smilarity matrix.csv'
    swscore_matrix = pd.read_csv(swpath, header=None).values

    Positives,Negatives = get_Positives_and_Negatives(positive_index, protein, RNA, swscore_matrix, savepath)


    Positives, Negatives_sort = get_edgelist(Positives, Negatives, 'sort', savepath,
                                             len(RNA)) 
    _, Negatives_random = get_edgelist(Positives, Negatives, 'random', savepath, len(RNA))
    _, Negatives_sort_random = get_edgelist(Positives, Negatives, 'sort_random', savepath, len(RNA))

    NPI_pos = np.zeros((len(RNA), len(protein)))
    NPI_pos[Positives.values[:, 0], Positives.values[:, 1]] = 1
    NPI_pos = pd.DataFrame(NPI_pos)
    NPI_pos.to_csv(savepath + 'NPI_pos.csv', index=False, header=None)

    NPI_neg_sort, NPI_neg_random, NPI_neg_sort_random = np.zeros((len(RNA), len(protein))), np.zeros(
        (len(RNA), len(protein))), np.zeros((len(RNA), len(protein)))
    NPI_neg_sort[Negatives_sort.values[:, 0], Negatives_sort.values[:, 1]] = 1
    NPI_neg_sort = pd.DataFrame(NPI_neg_sort)
    NPI_neg_sort.to_csv(savepath + 'NPI_neg_sort.csv', index=False, header=None)

    NPI_neg_random[Negatives_random.values[:, 0], Negatives_random.values[:, 1]] = 1
    NPI_neg_random = pd.DataFrame(NPI_neg_random)
    NPI_neg_random.to_csv(savepath + 'NPI_neg_random.csv', index=False, header=None)

    NPI_neg_sort_random[Negatives_sort_random.values[:, 0], Negatives_sort_random.values[:, 1]] = 1
    NPI_neg_sort_random = pd.DataFrame(NPI_neg_sort_random)
    NPI_neg_sort_random.to_csv(savepath + 'NPI_neg_sort_random.csv', index=False, header=None)

if __name__ == '__main__':


    get_NPInter('D:/CSU/RPI-MD-/RPI-MD/data/raw_data/NPInter10412_dataset.txt', 'D:/CSU/RPI-MD-/RPI-MD/data/generated_data/NPInter_10412/')

    get_RPI7317('D:/CSU/RPI-MD-/RPI-MD/data/raw_data/RPI7317.csv','D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI7317/')

    get_RPI2241('D:/CSU/RPI-MD-/RPI-MD/data/raw_data/','D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI2241/')

    get_NPInter4158('D:/CSU/RPI-MD-/RPI-MD/data/raw_data/','D:/CSU/RPI-MD-/RPI-MD/data/generated_data/NPInter_4158/')

    get_RPI369('D:/CSU/RPI-MD-/RPI-MD/data/raw_data/RPI369_all.txt','D:/CSU/RPI-MD-/RPI-MD/data/generated_data/RPI369/')



