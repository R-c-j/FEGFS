import numpy as np
import pandas as pd

def batch_read_GO_Term(GO_Term_path,expression_path,outpath):
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
    Deal_GO_Term_gene_T.reset_index(drop=True, inplace=True)
    Deal_GO_Term_gene_T = Deal_GO_Term_gene_T.T
    Deal_GO_Term_gene_T.reset_index(drop=True, inplace=True)
    Deal_GO_Term_gene_T = Deal_GO_Term_gene_T.T
    gene_ID_file = Deal_GO_Term_gene_T
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
    # GO_Term_ID_col = GO_Term_ID_col.append(0)
    GO_Term_gene = gene_ID_file.iloc[1:, :]

    Repetition_rate_sum_Frame = pd.DataFrame()
    for j in np.arange(len(GO_Term_ID)):
        Term1 = Term = gene_ID_file.iloc[1:, j].dropna()
        Term1 = list(Term1)
        Repetition_rate_sum = []
        # Start to expand according to one term, and double check with the remaining terms
        for i in np.arange(len(GO_Term_ID)):
            repeat_gene = []
            Term2 = gene_ID_file.iloc[1:, i].dropna()
            Term2 = list(Term2)
            repeat_gene = [x for x in Term1 if x in Term2]
            if len(Term1) >= len(Term2):
                Repetition_rate = len(repeat_gene) / len(Term2)
            else:
                Repetition_rate = len(repeat_gene) / len(Term1)
            Repetition_rate_sum.append(Repetition_rate)
        Repetition_rate_sum = np.array(Repetition_rate_sum)
        Repetition_rate_sum = pd.DataFrame(pd.DataFrame(Repetition_rate_sum).T)
        Repetition_rate_sum_Frame = pd.concat([Repetition_rate_sum_Frame, Repetition_rate_sum], axis=0)

    """The following operation is to generate a gene repetition rate matrix between terms"""
    Repetition_rate_sum_Frame = pd.DataFrame(Repetition_rate_sum_Frame)
    Repetition_rate_sum_Frame.reset_index(drop=True, inplace=True)
    Repetition_rate_sum_Frame1 = pd.concat([GO_Term_ID_row, Repetition_rate_sum_Frame], axis=1)
    Repetition_rate_sum_Frame1 = pd.DataFrame(Repetition_rate_sum_Frame1)
    Repetition_rate_sum_Frame1 = Repetition_rate_sum_Frame1.T
    Repetition_rate_sum_Frame1.reset_index(drop=True, inplace=True)
    Repetition_rate_sum_Frame1 = Repetition_rate_sum_Frame1.T
    Repetition_rate_sum_Frame2 = pd.concat([GO_Term_ID_col, Repetition_rate_sum_Frame1], axis=0)
    Repetition_rate_sum_Frame2 = pd.DataFrame(Repetition_rate_sum_Frame2)
    Repetition_rate_sum_Frame2.reset_index(drop=True, inplace=True)

    """The following operation is to change the repetition rate matrix into the upper triangular matrix"""
    Repetition_rate_matrix_np = np.array(Repetition_rate_sum_Frame2.iloc[1:, 1:])
    Repetition_rate_matrix_np = np.triu(Repetition_rate_matrix_np)
    # A diagonal matrix is generated to prepare for the subsequent matrix subtraction
    Diagonal_matrix = np.eye(Repetition_rate_matrix_np.shape[0])
    Repetition_rate_matrix_np = Repetition_rate_matrix_np - Diagonal_matrix
    Repetition_rate_matrix = pd.DataFrame(Repetition_rate_matrix_np)
    Repetition_rate_matrix_Frame1 = pd.concat([GO_Term_ID_row, Repetition_rate_matrix], axis=1)
    Repetition_rate_matrix_Frame1 = pd.DataFrame(Repetition_rate_matrix_Frame1)
    Repetition_rate_matrix_Frame1 = Repetition_rate_matrix_Frame1.T
    Repetition_rate_matrix_Frame1.reset_index(drop=True, inplace=True)
    Repetition_rate_matrix_Frame1 = Repetition_rate_matrix_Frame1.T
    Repetition_rate_matrix_Frame2 = pd.concat([GO_Term_ID_col, Repetition_rate_matrix_Frame1], axis=0)
    Repetition_rate_matrix_Frame2 = pd.DataFrame(Repetition_rate_matrix_Frame2)
    Repetition_rate_matrix_Frame2.reset_index(drop=True, inplace=True)
    Repetition_rate_matrix_Frame2 = pd.DataFrame(Repetition_rate_matrix_Frame2)

    Repetition_rate_matrix_col_row = Repetition_rate_matrix_Frame2.iloc[1:, 1:]
    row = list(Repetition_rate_matrix_Frame2.iloc[1:, 0])
    col = list(Repetition_rate_matrix_Frame2.iloc[0, 1:])
    Repetition_rate_than_50 = []
    row_col = []
    matrix_shape = Repetition_rate_matrix_col_row.shape[1]
    """The following operation is to screen out the term pairs with a repetition rate of 0.8"""
    for j in np.arange(matrix_shape):
        First_name_row = row[j]
        First = Repetition_rate_matrix_col_row.iloc[j, :]
        First = np.array(First)
        num = 0
        for i in np.arange(len(First)):
            if First[i] > 0.8:
                # Repetition_rate_than_50.append(First[i])
                First_name_col = col[i]
                row_col.append(First_name_row)
                row_col.append(First_name_col)
                # The repetition rate value is taken to two decimal places
                row_col.append(str(round(First[i], 2)))
                num += 1
        Repetition_rate_than_50.append(num)
    """The following operation is to generate a data table to record the repetition rate between pairs of terms"""
    Repetition_rate_than_50 = np.array(Repetition_rate_than_50)
    row_col = np.array(row_col)
    row_col = row_col.reshape(-1, 3)
    # Generate a data table of two term + repetition rate values
    row_col = pd.DataFrame(row_col)
    row_col.columns = pd.Series(['Term1', 'Term2', 'Repetition_rate'])
    row_col_1 = row_col.sort_values(by='Repetition_rate', ascending=False)
    row_col_1 = pd.DataFrame(row_col_1)
    row_col_rank = row_col_1.iloc[1:, :]
    row_col_rank.reset_index(drop=True, inplace=True)
    row_col_col3 = row_col_rank.iloc[0:, 2]
    row_col_col3 = np.array(row_col_col3)
    num = 0
    for i in np.arange(len(row_col_col3)):
        if row_col_col3[i] >= str(1):
            num += 1
    # Num1 function: count the number of repetition rate to 1, and extract the term pair with repetition rate of 1 in order to extract the term pair with repetition rate of 1 later
    row_col_rank_equal_1 = row_col_rank.iloc[0:num, :]

    row_col_equal_1_Term1 = row_col_rank_equal_1.iloc[:, 0]
    # row_col_equal_1_Term1 = pd.DataFrame(row_col_equal_1_Term1)
    row_col_equal_1_Term1.reset_index(drop=True, inplace=True)
    row_col_equal_1_Term2 = row_col_rank_equal_1.iloc[:, 1]
    # row_col_equal_1_Term2 = pd.DataFrame(row_col_equal_1_Term2)
    row_col_equal_1_Term2.reset_index(drop=True, inplace=True)

    Term_gene = gene_ID_file
    Term_name = Term_gene.iloc[0, :]
    Term_gene = Term_gene.iloc[1:, :]
    Term_gene.columns = np.array(Term_name)
    Term_gene.reset_index(drop=True, inplace=True)
    Term_gene_equal_1_summary = pd.DataFrame()
    Term_equal_1_name = []
    Term_equal_exclude_name = []
    for i in np.arange(len(row_col_equal_1_Term1)):
        # The gene names of term pairs were extracted
        Term1_equal_1_gene_1 = Term_gene[row_col_equal_1_Term1[i]].dropna()
        Term1_equal_1_gene_2 = pd.DataFrame(Term1_equal_1_gene_1)
        Term1_equal_1_gene = Term1_equal_1_gene_2.iloc[:, 0]
        Term2_equal_1_gene = Term_gene[row_col_equal_1_Term2[i]].dropna()
        Term2_equal_1_gene = pd.DataFrame(Term2_equal_1_gene)
        Term2_equal_1_gene = Term2_equal_1_gene.iloc[:, 0]
        if len(Term1_equal_1_gene) >= len(Term2_equal_1_gene):
            Term_equal_1_name.append(row_col_equal_1_Term1[i])
            Term_equal_exclude_name.append(row_col_equal_1_Term2[i])
        else:
            Term_equal_1_name.append(row_col_equal_1_Term2[i])
            Term_equal_exclude_name.append(row_col_equal_1_Term1[i])
        # Merge the genes in two terms
        Term_gene_equal_1_sum = pd.concat([Term1_equal_1_gene, Term2_equal_1_gene])
        # The repeated genes in the merged genes were removed
        Term_gene_equal_1_sum_1 = np.unique(np.array(Term_gene_equal_1_sum))
        Term_gene_equal_1_sum_2 = pd.DataFrame(Term_gene_equal_1_sum_1)
        # Add the merged terms to the new data table for summary
        Term_gene_equal_1_summary = pd.concat([Term_gene_equal_1_summary, Term_gene_equal_1_sum_2], axis=1)
    Term_gene_equal_1_summary = pd.DataFrame(Term_gene_equal_1_summary)

    Term_equal_1_name = pd.Series(Term_equal_1_name)
    # Term_equal_1_name = list(Term_equal_1_name)
    # Term_gene_equal_1_summary.columns = Term_equal_1_name
    # Term_equal_1_name = pd.DataFrame(Term_equal_1_name)
    Term_equal_1_name = pd.DataFrame(Term_equal_1_name).T
    Term_gene_equal_1 = Term_gene_equal_1_summary
    Term_gene_equal_1 = Term_gene_equal_1.T
    Term_gene_equal_1.reset_index(drop=True, inplace=True)
    Term_gene_equal_1 = Term_gene_equal_1.T
    Term_gene_equal_1 = pd.concat([Term_equal_1_name, Term_gene_equal_1], axis=0)

    Term_gene_equal_2 = pd.DataFrame(Term_gene_equal_1.T)
    # In order to extract term_ gene_ equal_ 2 to prepare for subsequent deletion of duplicate terms
    subset = list(Term_gene_equal_2)
    Term_gene_equal_2.drop_duplicates(subset=list(Term_gene_equal_2)[0], inplace=True)
    # Merge the terms with a repetition rate of 1 and delete the repeated terms
    Term_gene_equal_3 = pd.DataFrame(Term_gene_equal_2.T)

    GO_Term_LIST = list(GO_Term_ID)
    # row_col_equal_Term_LIST = list(Term_gene_equal_3.iloc[0,:])
    row_col_equal_1_Term_LIST = list(row_col_equal_1_Term1)
    row_col_equal_2_Term_LIST = list(row_col_equal_1_Term2)
    # Sort out the original term set and delete the term used in the merging from the original term set
    for i in range(len(GO_Term_LIST)):
        for j in GO_Term_LIST:
            if j in row_col_equal_1_Term_LIST:
                GO_Term_LIST.remove(j)
    for i in range(len(GO_Term_LIST)):
        for j in GO_Term_LIST:
            if j in row_col_equal_2_Term_LIST:
                GO_Term_LIST.remove(j)
    # After deleting from the original term name set, the remaining terms that do not participate in the merge
    GO_Term_LIST = pd.Series(GO_Term_LIST)
    GO_Term_LIST1 = Term_gene_equal_3.iloc[0, :]
    GO_Term_LIST_sum = pd.concat([GO_Term_LIST, GO_Term_LIST1])
    GO_Term_LIST_sum = pd.DataFrame(GO_Term_LIST_sum)
    gene_ID_file_T = gene_ID_file.T
    gene_ID_file_T_col_name = list(Term_gene_equal_2)
    # After merging the terms with a repetition rate of 1, delete the covered terms from the original term + gene set
    Term_gene_dealing = pd.merge(gene_ID_file_T, GO_Term_LIST_sum)
    Term_gene_dealing = Term_gene_dealing.T

    Term_gene_dealing = pd.DataFrame(Term_gene_dealing)
    Term_gene_dealing_1 = Term_gene_dealing.iloc[1:, :]
    Term_gene_dealing_1.columns = Term_gene_dealing.iloc[0, :]
    exclude_Term = []
    for i in np.arange(len(Term_equal_exclude_name)):
        if Term_equal_exclude_name[i] in list(Term_gene_dealing_1):
            exclude_Term.append(Term_equal_exclude_name[i])
    Term_gene_dealing_1.drop(columns=exclude_Term, axis=1, inplace=True)
    Term_gene_dealing_1.reset_index(drop=True, inplace=True)

    Term_gene_dealing_1_name_list = list(Term_gene_dealing_1)
    Term_gene_dealing_1_name = np.array(Term_gene_dealing_1_name_list)
    Term_gene_dealing_1_name = pd.DataFrame(Term_gene_dealing_1_name).T
    empty_list = []
    for i in np.arange(len(Term_gene_dealing_1_name_list)):
        empty_list.append(i)
    Term_gene_dealing_1.columns = empty_list
    Term_gene_dealing_1 = pd.concat([Term_gene_dealing_1_name, Term_gene_dealing_1], axis=0)
    Term_gene_dealing_1.reset_index(drop=True, inplace=True)

    gene_ID_file = Term_gene_dealing_1

    GO_Term_ID_row = pd.DataFrame(gene_ID_file.iloc[0, :])
    GO_Term_ID = np.array(gene_ID_file.iloc[0, :])
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
        Term1 = Term = gene_ID_file.iloc[1:, j].dropna()
        Term1 = list(Term1)
        Repetition_rate_sum = []

        for i in np.arange(len(GO_Term_ID)):
            Term2 = gene_ID_file.iloc[1:, i].dropna()
            Term2 = list(Term2)
            repeat_gene = [x for x in Term1 if x in Term2]
            if len(Term1) >= len(Term2):
                Repetition_rate = len(repeat_gene) / len(Term2)
            else:
                Repetition_rate = len(repeat_gene) / len(Term1)
            Repetition_rate_sum.append(Repetition_rate)
        Repetition_rate_sum = np.array(Repetition_rate_sum)
        Repetition_rate_sum = pd.DataFrame(pd.DataFrame(Repetition_rate_sum).T)
        Repetition_rate_sum_Frame = pd.concat([Repetition_rate_sum_Frame, Repetition_rate_sum], axis=0)

    Repetition_rate_sum_Frame = pd.DataFrame(Repetition_rate_sum_Frame)
    Repetition_rate_sum_Frame.reset_index(drop=True, inplace=True)
    Repetition_rate_sum_Frame1 = pd.concat([GO_Term_ID_row, Repetition_rate_sum_Frame], axis=1)
    Repetition_rate_sum_Frame1 = pd.DataFrame(Repetition_rate_sum_Frame1)
    Repetition_rate_sum_Frame1 = Repetition_rate_sum_Frame1.T
    Repetition_rate_sum_Frame1.reset_index(drop=True, inplace=True)
    Repetition_rate_sum_Frame1 = Repetition_rate_sum_Frame1.T
    Repetition_rate_sum_Frame2 = pd.concat([GO_Term_ID_col, Repetition_rate_sum_Frame1], axis=0)
    Repetition_rate_sum_Frame2 = pd.DataFrame(Repetition_rate_sum_Frame2)

    Repetition_rate_matrix_np = np.array(Repetition_rate_sum_Frame2.iloc[1:, 1:])
    Repetition_rate_matrix_np = np.triu(Repetition_rate_matrix_np)

    Diagonal_matrix = np.eye(Repetition_rate_matrix_np.shape[0])
    Repetition_rate_matrix_np = Repetition_rate_matrix_np - Diagonal_matrix
    Repetition_rate_matrix = pd.DataFrame(Repetition_rate_matrix_np)
    Repetition_rate_matrix_Frame1 = pd.concat([GO_Term_ID_row, Repetition_rate_matrix], axis=1)
    Repetition_rate_matrix_Frame1 = pd.DataFrame(Repetition_rate_matrix_Frame1)
    Repetition_rate_matrix_Frame1 = Repetition_rate_matrix_Frame1.T
    Repetition_rate_matrix_Frame1.reset_index(drop=True, inplace=True)
    Repetition_rate_matrix_Frame1 = Repetition_rate_matrix_Frame1.T
    Repetition_rate_matrix_Frame2 = pd.concat([GO_Term_ID_col, Repetition_rate_matrix_Frame1], axis=0)
    Repetition_rate_matrix_Frame2 = pd.DataFrame(Repetition_rate_matrix_Frame2)
    Repetition_rate_matrix_Frame2.reset_index(drop=True, inplace=True)
    Repetition_rate_matrix_Frame2 = pd.DataFrame(Repetition_rate_matrix_Frame2)

    Repetition_rate_matrix_col_row = Repetition_rate_matrix_Frame2.iloc[1:, 1:]
    row = list(Repetition_rate_matrix_Frame2.iloc[1:, 0])
    col = list(Repetition_rate_matrix_Frame2.iloc[0, 1:])
    Repetition_rate_than_50 = []
    row_col = []
    matrix_shape = Repetition_rate_matrix_col_row.shape[1]

    for j in np.arange(matrix_shape):
        First_name_row = row[j]
        First = Repetition_rate_matrix_col_row.iloc[j, :]
        First = np.array(First)
        num = 0
        for i in np.arange(len(First)):
            if float(First[i]) > float(0.8):
                # Repetition_rate_than_50.append(First[i])
                First_name_col = col[i]
                row_col.append(First_name_row)
                row_col.append(First_name_col)

                row_col.append(str(round(First[i], 2)))
                num += 1
        Repetition_rate_than_50.append(num)

    Repetition_rate_than_50 = np.array(Repetition_rate_than_50)
    row_col = np.array(row_col)
    row_col = row_col.reshape(-1, 3)

    row_col = pd.DataFrame(row_col)
    row_col.columns = pd.Series(['Term1', 'Term2', 'Repetition_rate'])
    row_col_1 = row_col.sort_values(by='Repetition_rate', ascending=False)

    row_col_equal_1_Term1 = row_col_1.iloc[:, 0]
    row_col_equal_1_Term2 = row_col_1.iloc[:, 1]

    Term_gene = gene_ID_file
    Term_name = Term_gene.iloc[0, :]
    Term_gene = Term_gene.iloc[1:, :]
    Term_gene.columns = np.array(Term_name)
    Term_gene.reset_index(drop=True, inplace=True)
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

            Term1_equal_1_gene_1 = Term_gene[row_col_equal_1_Term1[i]].dropna()
            Term1_equal_1_gene_2 = pd.DataFrame(Term1_equal_1_gene_1)
            Term1_equal_1_gene = Term1_equal_1_gene_2.iloc[:, 0]
            Term2_equal_1_gene = Term_gene[row_col_equal_1_Term2[i]].dropna()
            Term2_equal_1_gene = pd.DataFrame(Term2_equal_1_gene)
            Term2_equal_1_gene = Term2_equal_1_gene.iloc[:, 0]
            Term_gene_equal_1_sum = pd.concat([Term1_equal_1_gene, Term2_equal_1_gene])

            Term_gene_equal_1_sum_1 = np.unique(np.array(Term_gene_equal_1_sum))
            Term_gene_equal_1_sum_2 = pd.DataFrame(Term_gene_equal_1_sum_1)

            Term_gene_equal_1_summary = pd.concat([Term_gene_equal_1_summary, Term_gene_equal_1_sum_2], axis=1)
    Term_gene_equal_1_summary = pd.DataFrame(Term_gene_equal_1_summary)
    Term_equal_1_name = pd.Series(Term_equal_1_name)
    # Term_equal_1_name = pd.DataFrame(Term_equal_1_name).T
    Term_gene_equal_1_summary.columns = pd.Series(Term_equal_1_name)

    GO_Term_LIST = list(GO_Term_ID)
    row_col_equal_Term_LIST = repeat_Term

    for i in range(len(GO_Term_LIST)):
        for j in GO_Term_LIST:
            if j in row_col_equal_Term_LIST:
                GO_Term_LIST.remove(j)

    gene_ID_file = Term_gene_dealing_1
    gene_ID_file_name_1 = list(gene_ID_file)
    gene_ID_file_name = np.array(gene_ID_file_name_1)
    gene_ID_file_name = pd.DataFrame(gene_ID_file_name)
    gene_ID_file_name.reset_index(drop=True, inplace=True)
    gene_ID_file_name = gene_ID_file_name.T

    empty_list = []
    for i in np.arange(len(gene_ID_file_name_1)):
        empty_list.append(i)
    gene_ID_file.columns = empty_list
    gene_ID_file_1 = pd.concat([gene_ID_file_name, gene_ID_file], axis=0)
    gene_ID_file_1.reset_index(drop=True, inplace=True)
    GO_Term_extend = list(GO_Term_LIST) + list(Term_equal_1_name)
    GO_Term_extend = pd.Series(GO_Term_extend)
    GO_Term_LIST = pd.Series(GO_Term_LIST)
    GO_Term_LIST = pd.DataFrame(GO_Term_LIST)
    gene_ID_file_1_T = gene_ID_file_1.T
    gene_ID_file_1_T.reset_index(drop=True, inplace=True)
    gene_ID_file_1_T = gene_ID_file_1_T.iloc[:, 1:]
    gene_ID_file_1_T = gene_ID_file_1_T.T
    gene_ID_file_1_T.reset_index(drop=True, inplace=True)
    gene_ID_file_1_T = gene_ID_file_1_T.T
    Term_gene_dealing = pd.merge(gene_ID_file_1_T, GO_Term_LIST)
    Term_gene_dealing = Term_gene_dealing.T
    Term_gene_dealing_1 = Term_gene_dealing.iloc[1:, :]
    Term_gene_dealing_1.columns = list(Term_gene_dealing.iloc[0, :])
    Term_gene_equal_1_summary = pd.DataFrame(Term_gene_equal_1_summary)
    Term_gene_equal_1_summary.reset_index(drop=True, inplace=True)
    Term_gene_dealing_1 = pd.DataFrame(Term_gene_dealing_1)
    Term_gene_dealing_1.reset_index(drop=True, inplace=True)
    Term_gene_dealing_final = pd.concat([Term_gene_dealing_1, Term_gene_equal_1_summary], axis=1)

    GO_Term_gene_expression_path = expression_path
    outputpath = outpath
    gene_ID_file = Term_gene_dealing_final
    expression_file = pd.read_csv(GO_Term_gene_expression_path, encoding='UTF-8', low_memory=False)
    n = gene_ID_file.shape[1]
    for i in np.arange(n):
        # Represents the go class name for the following ID matching.
        ColNames = gene_ID_file.columns[i]
        # Remove the Na value in the go class. If not, an error will occur in the following vlookup matching operation.
        gene_ID_file_deal_NA = gene_ID_file[ColNames].dropna()
        gene_ID_file_deal_NA = pd.DataFrame(gene_ID_file_deal_NA)
        # The gene expression data matrix of go class was generated one by one by vlookup function matching operation.
        gene_ID_file_deal_NA_col = list(gene_ID_file_deal_NA)
        for j in np.arange(len(gene_ID_file_deal_NA_col)):
            gene_ID_file_deal_NA_1 = gene_ID_file_deal_NA.iloc[:, j]
            gene_ID_file_deal_NA_2 = pd.DataFrame(gene_ID_file_deal_NA_1)
            data_merge = pd.merge(left=expression_file, right=gene_ID_file_deal_NA_2, left_on="gene_ID",
                                  right_on=ColNames)
            # Remove the gene ID used in the last column matching after vlookup to ensure that only the gene ID of the first column is included.
            data_merge = data_merge.drop(ColNames, axis=1)
            data_merge = data_merge.drop('gene_ID', axis=1)
            data_merge.to_csv(outputpath + '/' + str(ColNames) + '_' + str(j) + '.csv', sep=',', index=False,
                              header=False)
    return


