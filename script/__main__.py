from batch_read_GO_Term import batch_read_GO_Term
from batch_read_GO_Term import batch_read_GO_Term_matrix
from Function_filter_Term_gene_repetitions_completely import screen_Term_gene_1
from Function_filter_Term_gene_repetitions import screen_Term_gene
from batch_single_exacting_feature_KPCA_method import batch_exacting_feature_KPCA
from clustering import clustering
import os
import time


""" Test """

def main_test():

    current_dir = os.path.dirname(os.path.abspath(__file__))

    GO_Term_path = current_dir + '/example/GO_Term.xlsx'
    outpath = current_dir + '/example/test_dealing.csv'
    batch_read_GO_Term(GO_Term_path=GO_Term_path, outputpath=outpath)

    Save_path = current_dir + '/example/test_duplicate_removal_1.csv'
    screen_Term_gene_1(GO_Term=outpath, Save_path=Save_path)

    Save_path_1 = current_dir + '/example/test_duplicate_removal.csv'
    screen_Term_gene(GO_Term=Save_path, Save_path=Save_path_1)

    GO_Term_gene_expression_path = current_dir + '/example/test_count_matrix.csv'
    outputpath1 = current_dir + '/example/Term_matrix'
    batch_read_GO_Term_matrix(GO_Term_gene_path=Save_path_1, GO_Term_gene_expression_path=GO_Term_gene_expression_path,
                              outputpath=outputpath1)

    outputpath2 = current_dir + '/example/Term_matrix'
    Save_path_feature_matrix = current_dir + '/example/test_feature_matrix.csv'
    batch_exacting_feature_KPCA(read_path=outputpath2, outputpath=Save_path_feature_matrix, n_components=0.4)

    label_path = current_dir + '/example/test_label.csv'
    clustering(read_path=Save_path_feature_matrix, n_clusters=5, label_path=label_path)


""" real data """

def main():
    # start = time.time()

    num = 5
    data = ['pollen_human', 'goolam', 'petal_human', 'biase', 'kelin', 'zeisel']
    cluster_num = [11, 5, 5, 4, 4, 7]
    n_component = [0.4, 0.4, 0.4, 0.4, 0.6, 0.8]
    data_name = data[num]
    n_components = n_component[num]
    n_clusters = cluster_num[num]

    GO_Term_path = 'GO_BP_MF_CC.xlsx'
    outputpath =  data_name + '_dealing.csv'
    batch_read_GO_Term(GO_Term_path, outputpath)

    GO_Term =  data_name + '_dealing.csv'
    Save_path = data_name + '_duplicate_removal_1.csv'
    screen_Term_gene_1(GO_Term, Save_path)

    GO_Term_1 =  data_name + '_duplicate_removal_1.csv'
    Save_path_1 =  data_name + '_duplicate_removal.csv'
    screen_Term_gene(GO_Term_1, Save_path_1)

    Save_path_2 =  data_name + '_duplicate_removal.csv'
    GO_Term_gene_expression_path =  data_name + '_count_matrix.csv'
    outputpath1 = 'Term_matrix'
    batch_read_GO_Term_matrix(GO_Term_gene_path=Save_path_2, GO_Term_gene_expression_path=GO_Term_gene_expression_path,
                              outputpath=outputpath1)

    read_path = 'Term_matrix'
    outputpath2 = data_name + '_KPCA_' + str(n_components) + '_cosine.csv'

    batch_exacting_feature_KPCA(read_path, outputpath2, n_components)

    read_path =  data_name + '_KPCA_' + str(n_components) + '_cosine.csv'
    label_path =  data_name + '_label.csv'
    clustering(n_clusters, read_path, label_path)

    # end = time.time()
    # print("+"*30)
    # print("Running time：%.4f second" %float(end - start))

if __name__ == '__main__':
    main()
    main_test()