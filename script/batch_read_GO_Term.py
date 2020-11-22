"""
This script file is divided into two parts of code:

The first part reads the genes contained in each term and generates the term gene expression matrix (CSV file)

In the second part, we match the generated term gene expression matrix with the original gene expression matrix, and generate a gene expression matrix (CSV file) for each term

"""
"""
Part I code
"""
import numpy as np
import pandas as pd
def batch_read_GO_Term(GO_Term_path,outputpath):
    """
    :param GO_Term_path:
    :param outputpath:
    :return:
    """
    GO_Term_path = GO_Term_path
    file = pd.read_excel(GO_Term_path, encoding='UTF-8')
    # Read and summarize terms_ ID
    GO_Term_sum_ID = file[['term_id']]
    # Read and summarize go_ Gene ID of term
    GO_Term_sum_gene = file[['intersections']]
    GO_Term_sum_gene_1 = pd.concat([GO_Term_sum_ID, GO_Term_sum_gene], axis=1)
    n = len(GO_Term_sum_gene_1)
    GO_Term_sum_gene_list = []

    # gene in GO term is added to the list
    for i in np.arange(n):
        GO_Term_sum_gene_list.append(GO_Term_sum_gene_1[i:i + 1])

    Deal_GO_Term_gene = pd.DataFrame()

    for i in np.arange(len(GO_Term_sum_gene_list)):

        Empty_Dataframe = pd.DataFrame()

        Dataframe_1 = np.array(GO_Term_sum_gene_list[i])
        Dataframe_2 = str(Dataframe_1)

        Dataframe_3 = Dataframe_2.split(' ')
        Dataframe_3 = str(Dataframe_3)
        Dataframe_3 = Dataframe_3.lstrip('["[[')
        Dataframe_3 = Dataframe_3.lstrip("'")
        Dataframe_3 = Dataframe_3.rstrip(']]"]')
        Dataframe_3 = Dataframe_3.rstrip("'")
        Dataframe_3 = Dataframe_3.split(',')
        Dataframe_4 = np.array(Dataframe_3)
        Dataframe_4_1 = str(Dataframe_4[0]).rstrip('\\n"')
        Dataframe_4_1 = Dataframe_4_1.rstrip("\'")
        Dataframe_4_1 = Dataframe_4_1.replace(':', "_")
        Dataframe_4[0] = Dataframe_4_1
        Dataframe_4_2 = str(Dataframe_4[2]).lstrip('"\ ')
        Dataframe_4_2 = Dataframe_4_2.lstrip("\'")
        Dataframe_4[2] = Dataframe_4_2
        for j in np.arange(len(Dataframe_4)):
            Empty_Dataframe[j] = [Dataframe_4[j]]

        Deal_GO_Term_gene = pd.concat([Deal_GO_Term_gene, Empty_Dataframe], axis=0)
        Deal_GO_Term_gene = Deal_GO_Term_gene.drop(columns=1)

        Deal_GO_Term_gene_T = pd.DataFrame(Deal_GO_Term_gene.T, index=Deal_GO_Term_gene.columns,
                                           columns=Deal_GO_Term_gene.index)
    Deal_GO_Term_gene_T.to_csv(str(outputpath), sep=',', index=False, header=False)
    return

def batch_read_GO_Term_matrix(GO_Term_gene_path,GO_Term_gene_expression_path,outputpath):

    gene_ID_file = pd.read_csv(GO_Term_gene_path, encoding='UTF-8',low_memory=False)
    expression_file = pd.read_csv(GO_Term_gene_expression_path, encoding='UTF-8',low_memory=False)
    n = gene_ID_file.shape[1]
    for i in np.arange(n):
        # Represents the go class name for the following ID matching.
        ColNames = gene_ID_file.columns[i]
        # Remove the Na value in the go class. If not, an error will occur in the following vlookup matching operation.
        gene_ID_file_deal_NA = gene_ID_file[ColNames].dropna()
        # The gene expression data matrix of go class was generated one by one by vlookup function matching operation.
        data_merge = pd.merge(left=expression_file, right=gene_ID_file_deal_NA, left_on="gene_ID", right_on=ColNames)
        # Remove the gene ID used in the last column matching after vlookup to ensure that only the gene ID of the first column is included.
        data_merge = data_merge.drop(ColNames, axis=1)
        data_merge = data_merge.drop('gene_ID', axis=1)
        data_merge.to_csv(outputpath+'/'+str(ColNames)+'.csv',sep=',', index=False,header=False)
    return



