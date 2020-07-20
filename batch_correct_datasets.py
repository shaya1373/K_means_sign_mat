import pandas as pd
import combat_master.combat as combat


def differentiate_same_col_names(input_df):
    """
    There might be multiple columns of the same name. Combat function raise an error in that case. This function makes sure
    that each column has a unique name.

    Args:
        input_df: pd DataFrame with repetitive column names
    Returns:
        pd DataFrame with unique column names, the same size as input_df
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


def remove_batch_effect(files_list, cell_type, address):
    """
    Args:
        files_list: List of strings: List of names of files that include cell_type
        cell_type: string:  Name of cell type we are looking for
    Returns:
        batch_corrected_data: (n_genes, n_samples) pd Dataframe:  Combined gene expression of cell types with removed
        batch effect
    """
    batch = pd.Series()
    data_dict = []
    for (i, file_name) in enumerate(files_list):
        data = pd.read_csv(address+file_name, sep="\t", header=0)
        data.set_index(list(data)[0], inplace=True)
        cell_type_data = data.loc[:, data.columns.str.contains(cell_type)]
        assert(cell_type_data.shape[1] != 0), cell_type + " is not present in " + file_name

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


def perform_batch_correction_on_all_files(cell_dict, address):
    batch_corrected_list = []
    for cell in cell_dict:
        file_list = cell_dict[cell]
        batch_corrected_list.append(remove_batch_effect(file_list, cell, address))

    batch_corrected_df = pd.concat(batch_corrected_list, axis=1, sort=True)
    batch_corrected_df.dropna(inplace=True)
    return batch_corrected_df


if __name__ == '__main__':
    address = "/Users/shaya/Desktop/University/Amherst/data_sets/Challehnge_data/Corrected/"
    cell_file_dict = {'CD8': ['ensembl_version_GSE98638.txt', 'ensembl_version_GSE114407.txt', 'GSE107011_Processed_data_TPM2.txt'],
                      'CD4': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE97861.txt', 'ensembl_version_GSE97862.txt',
                              'ensembl_version_GSE113891.txt', 'ensembl_version_GSE114407.txt', 'ensembl_version_GSE115978.txt'],
                      'B_': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE114407.txt'],
                      'mono': ['ensembl_version_GSE114407.txt', 'GSE107011_Processed_data_TPM2.txt'],
                      'NK': ['GSE107011_Processed_data_TPM2.txt', 'ensembl_version_GSE115978.txt'],
                      'Endo': ['ensembl_version_GSE102767.txt', 'ensembl_version_GSE113839.txt', 'ensembl_version_GSE115978.txt'],
                      'Fibro': ['ensembl_version_GSE113839.txt', 'GSE109448_rnaseq_gene_tpm.tsv', 'GSE109449_singlecell_rnaseq_gene_tpm.txt']}
    all_batch_corrected_df = perform_batch_correction_on_all_files(cell_file_dict, address)
    print(all_batch_corrected_df)
    neutro = pd.read_csv(address+"GSE107011_Processed_data_TPM2.txt", sep="\t", header=0)
    neutro.set_index(list(neutro)[0], inplace=True)
    neutro = neutro.loc[:, neutro.columns.str.contains('Neutro')]
    neutro = neutro.loc[~neutro.index.duplicated(keep='first')]

    data_out = pd.concat([all_batch_corrected_df, neutro], axis=1, sort=True)
    data_out.dropna(inplace=True)
    print(data_out)
    data_out.to_csv("batch_corrected_datasets.txt", sep="\t")
