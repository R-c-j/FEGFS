"""
Usage: t-SNE Dimension reduction
"""

from sklearn.manifold import TSNE
import pandas as pd
import numpy as np
import os

read_path = 'Term_matrix/'

os.chdir(read_path)

csv_name_list = os.listdir()

information_value = [2]
for j in np.arange(len(information_value)):

    TSNE_exact_feature_data_summary = pd.DataFrame()
    for i in np.arange(len(csv_name_list)):

        GO_gene_expression_data_path = 'Term_matrix/'+csv_name_list[i]

        GO_gene_expression_data = pd.read_csv(GO_gene_expression_data_path)
        GO_gene_expression_data_T = pd.DataFrame(GO_gene_expression_data.T, index=GO_gene_expression_data.columns,
                                           columns=GO_gene_expression_data.index)
        GO_gene_expression_data_T = GO_gene_expression_data_T.fillna(0)
        GO_gene_expression_data_ndarray = np.array(GO_gene_expression_data_T)

        TSNE_exact_feature_data = TSNE(n_components=information_value[j])

        TSNE_exact_feature_data = TSNE_exact_feature_data.fit_transform(GO_gene_expression_data_ndarray)
        TSNE_exact_feature_data = pd.DataFrame(TSNE_exact_feature_data)

        TSNE_exact_feature_data_summary=pd.concat([TSNE_exact_feature_data,TSNE_exact_feature_data_summary],axis=1)

    TSNE_exact_feature_data_summary_T = pd.DataFrame(TSNE_exact_feature_data_summary.T, index=TSNE_exact_feature_data_summary.columns,
                                           columns=TSNE_exact_feature_data_summary.index)

    outputpath = 'pollen_TSNE_'+str(information_value[j])+'.csv'

    TSNE_exact_feature_data_summary_T.to_csv(outputpath, sep=',', index=False, header=False)




