import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy import stats
import statsmodels.stats.multitest as multi


MAX_CLUSTERS = 10


def find_optimal_k(cell_expr, max_clusters):
    """
    Args:
        cell_expr: (n_samples, n_genes) pd DataFrame - Transposed DataFrame of cell expression
        max_clusters: integer - Specifies maximum number of of clusters to be considered
    Returns:
        k_opt: integer - optimal number of clusters based on the Silhouette method
    """
    sample_num = cell_expr.shape[0]
    k_vals = np.zeros(len(range(2, min(max_clusters, sample_num))))
    range_n_clusters = list(range(2, min(max_clusters, sample_num)))
    for i, n_clusters in enumerate(range_n_clusters):
        clusterer = KMeans(n_clusters=n_clusters, random_state=10, n_init=10, n_jobs=-1, tol=1e-3)
        cluster_labels = clusterer.fit_predict(cell_expr)
        silhouette_avg = silhouette_score(cell_expr, cluster_labels)
        k_vals[i] = silhouette_avg

    k_opt = np.argmax(k_vals)
    print(k_opt)
    return k_opt + 2


def cluster_each_cell(cells, cell_expr):
    """
    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - DataFrame of cell expression
    Returns:
        out_df: (n_genes, n_samples) pd DataFrame - clustered cell expression. each column shows the cell type and the
                cluster number that sample belongs to, e.g. CD4_subtype_1.
    """
    out_df = pd.DataFrame(index=cell_expr.index)
    for cell in cells:
        one_cell_expr = cell_expr.loc[:, cell_expr.columns.str.contains(cell)]

        assert (one_cell_expr.shape[1] != 0), cell + " is not present in the file"

        k_opt = find_optimal_k(one_cell_expr.T, MAX_CLUSTERS)
        clusterer = KMeans(n_clusters=k_opt, random_state=10, n_init=10, n_jobs=-1, tol=1e-3)
        cluster_labels = clusterer.fit_predict(one_cell_expr.T)

        one_cell_expr.columns = [cell + "_subtype_" + str(cluster_labels[i] + 1) + "." + str(i) for i in range(len(cluster_labels))]
        out_df[list(one_cell_expr)] = one_cell_expr
    out_df = out_df.rename_axis("genes")
    return out_df


def get_cluster_names(cells, cell_expr):
    """
    If the clustered dataset is imported from an external file, this function gets the number of clusters of
    each cell-type in the "cell" list, and assign name to each cluster in the format of cellName_subtype_clusterNum

    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - Transposed DataFrame of cell expression
    Returns:
        cluster_name_list: list of strings - list of all cells and clusters in the file in the format of CD4_subtype_1.
    """
    cluster_name_list = []
    for cell in cells:
        cell_type_list = list(cell_expr.loc[:, cell_expr.columns.str.contains(cell)])

        assert (len(cell_type_list) != 0), cell + " is not present in the file"

        cell_type_list = [cell_name[:cell_name.find(".")] for cell_name in cell_type_list]
        num_cell_cluster = len(np.unique(cell_type_list))
        for i in range(num_cell_cluster):
            cluster_name_list.append(cell + "_subtype_" + str(i+1))

    return cluster_name_list


def get_mean_of_each_cluster(cells, cell_expr):
    """
    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - Transposed DataFrame of cell expression
    Returns:
        mean_df: (n_genes, n_clusters) pd Dataframe - mean value of each cluster across all samples of that cluster.
    """
    cluster_name_list = get_cluster_names(cells, cell_expr)
    mean_df = pd.DataFrame(index=cell_expr.index)
    for cluster_name in cluster_name_list:
        one_cell_mean = all_cells.loc[:, all_cells.columns.str.contains(cluster_name)].mean(axis=1)
        mean_df[cluster_name] = one_cell_mean
    return mean_df


def detect_highly_expr_genes(one_cluster, other_clusters):
    """
    Args:
        one_cluster: (n_genes, n_samples) pd DataFrame- A single cluster with its samples
        other_clusters: (n_genes, n_samples) pd DataFrame - Other clusters!
    Returns:
        list of strings - returns the first 100 genes that are highly expressed in one_cluster, but lowly expressed in
        other_clusters
    """
    one_cluster = np.log2(one_cluster + 1)
    other_clusters = np.log2(other_clusters + 1)
    one_cluster_min = one_cluster.min(axis=1)
    other_clusters_mean = other_clusters.mean(axis=1)
    diff = one_cluster_min - other_clusters_mean
    sub_df = diff[diff > 1]

    sub_df = sub_df.sort_values(ascending=False)
    if len(sub_df) > 100:
        sub_df = sub_df.iloc[:100]
    return list(sub_df.index)


def get_differentially_expr_genes(cells, cell_expr, mean_df, p_val_thresh=0.05):
    """
    Args:
        cells: list of strings - Specifies what type of cells we have
        cell_expr: (n_genes, n_samples) pd DataFrame - Transposed DataFrame of cell expression
        mean_df: (n_genes, n_clusters) pd Dataframe - mean value of each cluster across all samples of that cluster.
        p_val_thresh: float - Threshold for adjuster p value
    Returns:
        signature matrix - pd DataFrame
    """
    cluster_name_list = get_cluster_names(cells, cell_expr)
    final_selected_genes = []
    for cluster_name in cluster_name_list:
        one_cell_expr = all_cells.loc[:, all_cells.columns.str.contains(cluster_name)]
        other_cell_exp = all_cells.loc[:, ~all_cells.columns.str.contains(cluster_name)]

        _, p_val = stats.ttest_ind(one_cell_expr, other_cell_exp, axis=1, equal_var=False)
        _, q_val, _, _ = multi.multipletests(p_val)
        one_cell_expr = one_cell_expr[q_val < p_val_thresh]
        other_cell_exp = other_cell_exp[q_val < p_val_thresh]
        selected_genes = detect_highly_expr_genes(one_cell_expr, other_cell_exp)
        final_selected_genes += selected_genes
    final_selected_genes = np.unique(final_selected_genes)

    return mean_df.loc[final_selected_genes]


if __name__ == "__main__":
    cluster_flag = False
    p_val_flag = True
    cell_types = ["CD8", "CD4", "B_cell", "NK", "mono", "Endothelial", "Fibroblast", "Neutrophils"]

    if cluster_flag:
        all_cells = pd.read_csv("batch_corrected_datasets.txt", sep="\t")
        all_cells.set_index(list(all_cells)[0], inplace=True)
        clustered_cells = cluster_each_cell(cell_types, all_cells)
        clustered_cells.to_csv("batch_corrected_datasets_clustered.txt", sep="\t")

    if p_val_flag:
        all_cells = pd.read_csv("batch_corrected_datasets_clustered.txt", sep="\t")
        all_cells.set_index(list(all_cells)[0], inplace=True)
        all_cells[all_cells < 0] = 0

        mean_of_cluster = get_mean_of_each_cluster(cell_types, all_cells)
        final_sign_mat = get_differentially_expr_genes(cell_types, all_cells, mean_of_cluster)
        final_sign_mat.to_csv("K_means_signature_matrix_q_val.txt", sep="\t")
