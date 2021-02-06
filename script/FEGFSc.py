from batch_read_GO_Term import batch_read_GO_Term
from batch_read_GO_Term import batch_read_GO_Term_matrix
from Function_filter_Term_gene_repetitions_completely import screen_Term_gene_1
from Function_filter_Term_gene_repetitions import screen_Term_gene
from batch_single_exacting_feature_KPCA_method import batch_exacting_feature_KPCA
from clustering import clustering
import argparse
import os
import time

parser = argparse.ArgumentParser(description="Implementation of FEGFS")
parser.add_argument('-d','--data_name',dest='data_name',type=str)
parser.add_argument('-c','--n_clusters',dest='n_clusters',type=int)
parser.add_argument('-n','--n_components',dest='n_components',type=float,default=0.4)
parser.add_argument('-i','--input_GO_Term_path',dest='GO_Term_path',type=str)
parser.add_argument('-e','--expression_path',dest='expression_path',type=str)
parser.add_argument('-o','--outputpath1',dest='outputpath1',type=str)
parser.add_argument('-l','--label_path',dest='label_path',type=str)
args = parser.parse_args()

data_name = args.data_name
n_components = args.n_components
n_clusters = args.n_clusters

GO_Term_path = args.GO_Term_path
outputpath =  data_name + '_dealing.csv'
batch_read_GO_Term(GO_Term_path, outputpath)

GO_Term =  outputpath
Save_path = data_name + '_duplicate_removal_1.csv'
screen_Term_gene_1(GO_Term, Save_path)

GO_Term_1 =  Save_path
Save_path_1 =  data_name + '_duplicate_removal.csv'
screen_Term_gene(GO_Term_1, Save_path_1)

Save_path_2 =  Save_path_1
GO_Term_gene_expression_path = args.expression_path
outputpath1 = args.outputpath1
batch_read_GO_Term_matrix(GO_Term_gene_path=Save_path_2, GO_Term_gene_expression_path=GO_Term_gene_expression_path,
                          outputpath=outputpath1)

read_path = outputpath1
outputpath2 = data_name + '_KPCA_' + str(n_components) + '_cosine.csv'

batch_exacting_feature_KPCA(read_path, outputpath2, n_components)

read_path =  outputpath2
label_path =  args.label_path
clustering(n_clusters, read_path, label_path)

