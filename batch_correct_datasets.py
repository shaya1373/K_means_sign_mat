import pandas as pd
import numpy as np
from sklearn.svm import NuSVR
import combat_master.combat as combat


def differentiate_same_col_names(input_df):
    """
    param input_df: DataFrame with repetitive column names
    return: DataFrame with unique column names
    """
    col_names = list(input_df)
    col_list = []
    col_dict = {}
    for col_name in col_names:
        if col_name in col_dict:
            col_dict[col_name] += 1
            col_list.append(col_name + "_" + str(col_dict[col_name]))
        else:
            col_dict[col_name] = 0
            col_list.append(col_name)
    input_df.columns = col_list
    return input_df


def remove_batch_effect(files_dict, cell_type):
    """
    param files_dict: List of names of files that include cell_type (list)
    param cell_type: Name of cell type we are looking for (str)
    return: Combined gene expression of cell types with removed batch effect (DataFrame: num_genes * num_samples)
    """
    batch = pd.Series()
    data_dict = []
    address = "/Users/shaya/Desktop/University/Amherst/data_sets/Challehnge_data/Corrected/"
    for (i, file_name) in enumerate(files_dict):
        data = pd.read_csv(address+file_name, sep="\t", header=0)
        data.set_index(list(data)[0], inplace=True)
        cell_type_data = data.loc[:, data.columns.str.contains(cell_type)]
        temp_batch = pd.Series([i for _ in range(len(list(cell_type_data)))], index=list(cell_type_data))
        batch = batch.append(temp_batch)
        cell_type_data = cell_type_data.loc[~cell_type_data.index.duplicated(keep='first')]
        data_dict.append(cell_type_data)

    data_combined = pd.concat(data_dict, axis=1, sort=True)
    data_combined = data_combined.apply(lambda row: row.fillna(row.mean()), axis=1)
    data_combined.dropna(inplace=True)
    var = data_combined.var(axis=1)
    data_combined = data_combined[var > 0]
    data_combined = differentiate_same_col_names(data_combined)
    batch.index = list(data_combined)
    batch_corrected_data = combat.combat(data_combined, batch)
    return batch_corrected_data


if __name__ == '__main__':
    cell_file_dict = {'CD8': ['ensembl_version_GSE98638.txt', 'ensembl_version_GSE114407.txt', 'GSE107011_Processed_data_TPM2.txt'],
                      'CD4': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE97861.txt', 'ensembl_version_GSE97862.txt',
                              'ensembl_version_GSE113891.txt', 'ensembl_version_GSE114407.txt', 'ensembl_version_GSE115978.txt'],
                      'B_cell': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE114407.txt'],
                      'mono': ['ensembl_version_GSE114407.txt', 'GSE107011_Processed_data_TPM2.txt'],
                      'NK': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE115978.txt'],
                      'Endo': ['ensembl_version_GSE102767.txt', 'ensembl_version_GSE113839.txt', 'ensembl_version_GSE115978.txt'],
                      'Fibro': ['ensembl_version_GSE113839.txt', 'GSE109448_rnaseq_gene_tpm.tsv', 'GSE109449_singlecell_rnaseq_gene_tpm.txt']}

    CD8_batch_corrected = remove_batch_effect(["ensembl_version_GSE98638.txt", "ensembl_version_GSE114407.txt",
                         "GSE107011_Processed_data_TPM2.txt"], "CD8")
    CD4_batch_corrected = remove_batch_effect(["GSE107011_Processed_data_TPM2.txt", "ensembl_version_GSE97861.txt",
                               "ensembl_version_GSE97862.txt", "ensembl_version_GSE113891.txt",
                               "ensembl_version_GSE114407.txt", "ensembl_version_GSE115978.txt"], "CD4")
    B_batch_corrected = remove_batch_effect(["GSE107011_Processed_data_TPM2.txt", "ensembl_version_GSE114407.txt",
                               "ensembl_version_GSE115978.txt"], "B_")
    mono_batch_corrected = remove_batch_effect(["ensembl_version_GSE114407.txt", "GSE107011_Processed_data_TPM2.txt"],
                                               "mono")
    NK_batch_corrected = remove_batch_effect(["GSE107011_Processed_data_TPM2.txt", "ensembl_version_GSE115978.txt"],
                                               'NK')
    endo_batch_corrected = remove_batch_effect(["ensembl_version_GSE102767.txt", "ensembl_version_GSE113839.txt",
                               "ensembl_version_GSE115978.txt"], "Endo")
    fibro_batch_corrected = remove_batch_effect(["ensembl_version_GSE113839.txt", "GSE109448_rnaseq_gene_tpm.tsv",
                               "GSE109449_singlecell_rnaseq_gene_tpm.txt"], "Fibro")
    address1 = "/Users/shaya/Desktop/University/Amherst/data_sets/Challehnge_data/Corrected/"
    neutro = pd.read_csv(address1+"GSE107011_Processed_data_TPM2.txt", sep="\t", header=0)
    neutro.set_index(list(neutro)[0], inplace=True)
    neutro = neutro.loc[:, neutro.columns.str.contains('Neutro')]
    neutro = neutro.loc[~neutro.index.duplicated(keep='first')]

    data_out = pd.concat([CD4_batch_corrected, CD8_batch_corrected, B_batch_corrected, mono_batch_corrected,
                          NK_batch_corrected, endo_batch_corrected, fibro_batch_corrected, neutro], axis=1, sort=True)
    data_out.dropna(inplace=True)
    print(data_out)
    data_out.to_csv("batch_corrected_datasets.txt", sep="\t")
