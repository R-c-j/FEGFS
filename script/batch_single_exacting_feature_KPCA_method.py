
from sklearn.decomposition import KernelPCA
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import os
import time

def batch_exacting_feature_KPCA(read_path,outputpath,n_components):
    """
    This method is suitable for reading and dimensionality reduction of term gene expression data sets after go enrichment analysis.
    :param read_path: Read file path
    :param outputpath: Save path of characteristic matrix (CSV file)
    :param n_components: The proportion of principal components to be retained in feature extraction [0.8,0.95]
    :param kernel: Kernel method, example:'linear','poly','rbf','sigmoid','cosine'
    :param gamma: Parameter gamma value: int type, the general best is: 15
    :return: PCA characteristic data matrix
    """
    os.chdir(read_path)
    csv_name_list = os.listdir()
    information_value = [n_components]
    for j in np.arange(len(information_value)):

        KPCA_exact_feature_data_summary = pd.DataFrame()
        for i in np.arange(len(csv_name_list)):
            GO_gene_expression_data_path = read_path +'/'+csv_name_list[i]
            GO_gene_expression_data = pd.read_csv(GO_gene_expression_data_path, header=None)
            GO_gene_expression_data_T = pd.DataFrame(GO_gene_expression_data.T, index=GO_gene_expression_data.columns,
                                                     columns=GO_gene_expression_data.index)
            GO_gene_expression_data_T = GO_gene_expression_data_T.fillna('0')
            GO_gene_expression_data_ndarray = np.array(GO_gene_expression_data_T)

            pca = PCA(n_components=information_value[j])
            PCA_exact_feature_data = pca.fit_transform(GO_gene_expression_data_ndarray)
            value_PCA_exact_feature_data = PCA_exact_feature_data.shape[1]
            KPCA_exact_feature_data = KernelPCA(n_components=value_PCA_exact_feature_data, kernel='cosine')
            KPCA_exact_feature_data = KPCA_exact_feature_data.fit_transform(GO_gene_expression_data_ndarray)
            KPCA_exact_feature_data = pd.DataFrame(KPCA_exact_feature_data)
            KPCA_exact_feature_data_summary = pd.concat([KPCA_exact_feature_data, KPCA_exact_feature_data_summary],
                                                        axis=1)
        KPCA_exact_feature_data_summary_T = pd.DataFrame(KPCA_exact_feature_data_summary.T,
                                                         index=KPCA_exact_feature_data_summary.columns,
                                                         columns=KPCA_exact_feature_data_summary.index)
        KPCA_exact_feature_data_summary_T.to_csv(str(outputpath), sep=',', index=False, header=False)


def single_exacting_feature_KPCA(read_path,outputpath,n_components,kernel):
    """
    usage:It is suitable for reading gene expression matrix directly and reducing dimension
    :param read_path: Read file path
    :param outputpath: Save path of characteristic matrix (CSV file)
    :param n_components: The proportion of principal components to be retained in feature extraction [0.8,0.95]
    :param kernel: Kernel method, example:'linear','poly','rbf','sigmoid','cosine'
    :param gamma: Parameter gamma value: int type, the general best is: 15
    :return: PCA characteristic data matrix
    """

    GO_gene_expression_data_path = read_path
    GO_gene_expression_data = pd.read_csv(GO_gene_expression_data_path,header=0,index_col=0)
    GO_gene_expression_data = GO_gene_expression_data.fillna(0)
    GO_gene_expression_data_T = pd.DataFrame(GO_gene_expression_data.T, index=GO_gene_expression_data.columns,
                                             columns=GO_gene_expression_data.index)
    GO_gene_expression_data_ndarray = np.array(GO_gene_expression_data_T)
    pca = PCA(n_components=n_components)
    PCA_exact_feature_data = pca.fit_transform(GO_gene_expression_data_ndarray)
    value_PCA_exact_feature_data = PCA_exact_feature_data.shape[1]
    KPCA_exact_feature_data = KernelPCA(n_components=value_PCA_exact_feature_data,kernel=kernel)
    KPCA_exact_feature_data = KPCA_exact_feature_data.fit_transform(GO_gene_expression_data_ndarray)
    KPCA_exact_feature_data = pd.DataFrame(KPCA_exact_feature_data)
    KPCA_exact_feature_data_T = pd.DataFrame(KPCA_exact_feature_data.T,
                                                index=KPCA_exact_feature_data.columns,
                                                columns=KPCA_exact_feature_data.index)
    KPCA_exact_feature_data_T.to_csv(outputpath, sep=',', index=False, header=False)


