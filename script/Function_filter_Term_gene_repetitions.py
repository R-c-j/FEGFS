import pandas as pd
import numpy as np

def screen_Term_gene(GO_Term,Save_path):
    """
    :param data_name: pollen
    :param GO_Term: test_duplicate_removal_1.csv
    :param Save_path: test_duplicate_removal.csv
    :param repeat_value: 0.9
    :return
    """
    gene_ID_file = pd.read_csv(str(GO_Term), header=None)

    GO_Term_ID = gene_ID_file.iloc[0, :]
    GO_Term_ID_row = pd.DataFrame(GO_Term_ID)
    GO_Term_ID_col = pd.DataFrame(GO_Term_ID_row.T)
    data = np.arange(1).reshape(1, 1)
    Dataframe = pd.DataFrame(data)
    Dataframe[0] = None
    GO_Term_ID_col = pd.concat([Dataframe, GO_Term_ID_col], axis=1)
    GO_Term_ID_col = GO_Term_ID_col.T
    GO_Term_ID_col.reset_index(drop=True, inplace=True)
    GO_Term_ID_col = GO_Term_ID_col.T

    GO_Term_gene = gene_ID_file.iloc[1:, :]

    Repetition_rate_sum_Frame = pd.DataFrame()
    for j in np.arange(len(GO_Term_ID)):
        Term1 = Term = gene_ID_file.iloc[1:,j].dropna()
        Term1 = list(Term1)
        Repetition_rate_sum = []
        for i in np.arange(len(GO_Term_ID)):
            Term2 = gene_ID_file.iloc[1:,i].dropna()
            Term2 = list(Term2)
            repeat_gene = [x for x in Term1 if x in Term2]
            if len(Term1)>=len(Term2):
                Repetition_rate = len(repeat_gene)/len(Term2)
            else:
                Repetition_rate = len(repeat_gene) / len(Term1)
            Repetition_rate_sum.append(Repetition_rate)
        Repetition_rate_sum = np.array(Repetition_rate_sum)
        Repetition_rate_sum = pd.DataFrame(pd.DataFrame(Repetition_rate_sum).T)
        Repetition_rate_sum_Frame = pd.concat([Repetition_rate_sum_Frame,Repetition_rate_sum],axis=0)

    """The following operation is to generate a gene repetition rate matrix between terms"""

    Repetition_rate_sum_Frame = pd.DataFrame(Repetition_rate_sum_Frame)
    Repetition_rate_sum_Frame.reset_index(drop=True, inplace=True)
    Repetition_rate_sum_Frame1 = pd.concat([GO_Term_ID_row,Repetition_rate_sum_Frame],axis=1)
    Repetition_rate_sum_Frame1 = pd.DataFrame(Repetition_rate_sum_Frame1)
    Repetition_rate_sum_Frame1 = Repetition_rate_sum_Frame1.T
    Repetition_rate_sum_Frame1.reset_index(drop=True, inplace=True)
    Repetition_rate_sum_Frame1 = Repetition_rate_sum_Frame1.T
    Repetition_rate_sum_Frame2 = pd.concat([GO_Term_ID_col,Repetition_rate_sum_Frame1],axis=0)
    Repetition_rate_sum_Frame2 = pd.DataFrame(Repetition_rate_sum_Frame2)

    """The following operation is to change the repetition rate matrix into the upper triangular matrix"""

    Repetition_rate_matrix_np = np.array(Repetition_rate_sum_Frame2.iloc[1:,1:])
    Repetition_rate_matrix_np = np.triu(Repetition_rate_matrix_np)
    #A diagonal matrix is generated to prepare for the subsequent matrix subtraction
    Diagonal_matrix = np.eye(Repetition_rate_matrix_np.shape[0])
    Repetition_rate_matrix_np = Repetition_rate_matrix_np-Diagonal_matrix
    Repetition_rate_matrix = pd.DataFrame(Repetition_rate_matrix_np)
    Repetition_rate_matrix_Frame1 = pd.concat([GO_Term_ID_row,Repetition_rate_matrix],axis=1)
    Repetition_rate_matrix_Frame1 = pd.DataFrame(Repetition_rate_matrix_Frame1)
    Repetition_rate_matrix_Frame1 = Repetition_rate_matrix_Frame1.T
    Repetition_rate_matrix_Frame1.reset_index(drop=True, inplace=True)
    Repetition_rate_matrix_Frame1 = Repetition_rate_matrix_Frame1.T
    Repetition_rate_matrix_Frame2 = pd.concat([GO_Term_ID_col,Repetition_rate_matrix_Frame1],axis=0)
    Repetition_rate_matrix_Frame2 = pd.DataFrame(Repetition_rate_matrix_Frame2)
    Repetition_rate_matrix_Frame2.reset_index(drop=True, inplace=True)
    Repetition_rate_matrix_Frame2 = pd.DataFrame(Repetition_rate_matrix_Frame2)

    # Get row name, column name
    Repetition_rate_matrix_col_row = Repetition_rate_matrix_Frame2.iloc[1:, 1:]
    row = list(Repetition_rate_matrix_Frame2.iloc[1:, 0])
    col = list(Repetition_rate_matrix_Frame2.iloc[0, 1:])
    Repetition_rate_than_50 = []
    row_col = []
    matrix_shape = Repetition_rate_matrix_col_row.shape[1]
    """The following operation is to screen out the term pairs with a repetition rate of 0.9"""
    for j in np.arange(matrix_shape):
        First_name_row = row[j]
        First = Repetition_rate_matrix_col_row.iloc[j, :]
        First = np.array(First)
        num = 0
        for i in np.arange(len(First)):
            if First[i] > float(0.8):
                # Repetition_rate_than_50.append(First[i])
                First_name_col = col[i]
                row_col.append(First_name_row)
                row_col.append(First_name_col)

                row_col.append(str(round(First[i], 2)))
                num += 1
        Repetition_rate_than_50.append(num)

    """The following operation is to generate a data table to record the repetition rate between pairs of terms"""
    Repetition_rate_than_50 = np.array(Repetition_rate_than_50)
    row_col = np.array(row_col)
    row_col = row_col.reshape(-1, 3)

    row_col = pd.DataFrame(row_col)
    row_col.columns = pd.Series(['Term1', 'Term2', 'Repetition_rate'])
    row_col_1 = row_col.sort_values(by='Repetition_rate', ascending=False)

    row_col_equal_1_Term1 = row_col_1.iloc[:, 0]
    row_col_equal_1_Term2 = row_col_1.iloc[:, 1]

    Term_gene = pd.read_csv(str(GO_Term), header=0)
    Term_equal_1_name = []
    Term_gene_equal_1_summary = pd.DataFrame()
    repeat_Term = []
    for i in np.arange(len(row_col_equal_1_Term1)):
        if (row_col_equal_1_Term1[i] in repeat_Term) or (row_col_equal_1_Term2[i] in repeat_Term):
            continue
        else:
            Term_equal_1_name.append(row_col_equal_1_Term1[i] + '+' + row_col_equal_1_Term2[i])
            repeat_Term.append(row_col_equal_1_Term1[i])
            repeat_Term.append(row_col_equal_1_Term2[i])
            # The gene names of term pairs were extracted
            Term1_equal_1_gene = Term_gene[row_col_equal_1_Term1[i]].dropna()
            Term2_equal_1_gene = Term_gene[row_col_equal_1_Term2[i]].dropna()
            Term_gene_equal_1_sum = pd.concat([Term1_equal_1_gene, Term2_equal_1_gene])
            # The repeated genes in the merged genes were removed
            Term_gene_equal_1_sum_1 = np.unique(np.array(Term_gene_equal_1_sum))
            Term_gene_equal_1_sum_2 = pd.DataFrame(Term_gene_equal_1_sum_1)
            # Add the merged terms to the new data table for summary
            Term_gene_equal_1_summary = pd.concat([Term_gene_equal_1_summary, Term_gene_equal_1_sum_2], axis=1)
    Term_gene_equal_1_summary = pd.DataFrame(Term_gene_equal_1_summary)
    Term_equal_1_name = pd.Series(Term_equal_1_name)

    Term_gene_equal_1_summary.columns = pd.Series(Term_equal_1_name)
    """The following operation is to delete the terms with a repetition rate of more than 0.9 from the original term name set"""
    GO_Term_LIST = list(GO_Term_ID)
    row_col_equal_Term_LIST = repeat_Term
    for i in range(len(GO_Term_LIST)):
        for j in GO_Term_LIST:
            if j in row_col_equal_Term_LIST:
                GO_Term_LIST.remove(j)

    GO_Term_extend = GO_Term_LIST + list(Term_equal_1_name)
    GO_Term_extend = pd.Series(GO_Term_extend)
    GO_Term_LIST = pd.Series(GO_Term_LIST)
    GO_Term_LIST = pd.DataFrame(GO_Term_LIST)
    gene_ID_file_T = gene_ID_file.T
    Term_gene_dealing = pd.merge(gene_ID_file_T, GO_Term_LIST)
    Term_gene_dealing = Term_gene_dealing.T
    Term_gene_dealing_1 = Term_gene_dealing.iloc[1:, :]
    Term_gene_dealing_1.columns = list(Term_gene_dealing.iloc[0, :])
    Term_gene_equal_1_summary = pd.DataFrame(Term_gene_equal_1_summary)
    Term_gene_equal_1_summary.reset_index(drop=True, inplace=True)
    Term_gene_dealing_1 = pd.DataFrame(Term_gene_dealing_1)
    Term_gene_dealing_1.reset_index(drop=True, inplace=True)
    Term_gene_dealing_final = pd.concat([Term_gene_dealing_1, Term_gene_equal_1_summary], axis=1)
    return Term_gene_dealing_final.to_csv(str(Save_path), sep=',',index=False, header=True)









































