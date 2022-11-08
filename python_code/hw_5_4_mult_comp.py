import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest


def get_names(table):
    return table.columns[:-1].values


def check_intervals_intersect(first_ci: tuple, second_ci: tuple):
    are_intersect = False
    if (first_ci[1] > second_ci[0]) and (first_ci[0] < second_ci[1]):
        are_intersect = True
    return are_intersect


def interval_95(table, length):
    ans = st.t.interval(confidence=0.95,
                        df=length - 1,
                        loc=np.mean(table),
                        scale=st.sem(table))
    return ans


def interval_to_tables(first_table, second_table):
    first_intervals = []
    second_intervals = []
    gene_names = get_names(first_table)
    for gene in gene_names:
        df_f_length = len(first_table[gene])
        first_int = interval_95(first_table[gene], df_f_length)
        first_intervals.append(first_int)
        df_s_length = len(second_table[gene])
        second_int = interval_95(second_table[gene], df_s_length)
        second_intervals.append(second_int)
    return first_intervals, second_intervals


def check_dge_with_ci(first_table, second_table):
    f_intervals, s_intervals = interval_to_tables(first_table, second_table)
    ci_test_results = []
    for f_int, s_int in zip(
            f_intervals, s_intervals
    ):
        ci_test_results.append(check_intervals_intersect(f_int, s_int))
    return ci_test_results


def corr_mult_comp(p_values: list, method: str) -> list:
    from statsmodels.stats.multitest import multipletests as mpt
    return list(mpt(p_values, method=method)[1])


def find_p_val(first_table, second_table):
    gene_names = get_names(first_table)
    z_stat_val = []
    for gene in gene_names:
        ans = ztest(
            first_table[gene],
            second_table[gene]
        )
        z_stat_val.append(ans)
    z_test_p_values = []
    for num in z_stat_val:
        z_test_p_values.append(num[1])
    return z_test_p_values


def test_p_val(p_values):
    z_test_results = []
    for val in p_values:
        if val < 0.05:
            z_test_results.append(True)
        else:
            z_test_results.append(False)
    return z_test_results


def check_dge_with_ztest(first_table, second_table, results=True, values=False, method=None):
    z_test_p_values = find_p_val(first_table, second_table)
    z_test_results = test_p_val(z_test_p_values)
    if method:
        cor_p_values = corr_mult_comp(z_test_p_values, method=method)
        z_test_cor_res = test_p_val(cor_p_values)
        if results:
            return z_test_cor_res
        if values:
            return cor_p_values
    else:
        if results:
            return z_test_results
        if values:
            return z_test_p_values


def demonstrate_clt(expressions, sample_size: int = 1000, n_samples: int = 1000):
    randoms = np.array([
        np.random.choice(
            expressions, size=sample_size
        ) for _ in range(n_samples)
    ])
    mean_expressions = randoms.mean()
    return mean_expressions


def mean_differ(first_table, second_table):
    gene_names = get_names(first_table)
    f_mean = list(
        map(
            lambda gene: demonstrate_clt(first_table[gene], sample_size=100, n_samples=100),
            gene_names,
        )
    )
    s_mean = list(
        map(
            lambda gene: demonstrate_clt(second_table[gene], sample_size=100, n_samples=100),
            gene_names,
        )
    )
    mean_diff = np.subtract(f_mean, s_mean)
    return mean_diff


def create_table(first_cell_type_expressions_path, second_cell_type_expressions_path, save_results_table, method=None):
    first_table = pd.read_csv(f"{first_cell_type_expressions_path}", index_col=0)
    second_table = pd.read_csv(f"{second_cell_type_expressions_path}", index_col=0)
    ci_test_results = check_dge_with_ci(first_table, second_table)
    z_test_results = check_dge_with_ztest(first_table, second_table, results=True)
    z_test_p_values = check_dge_with_ztest(first_table, second_table, results=False, values=True)
    mean_diff = mean_differ(first_table, second_table)
    results = {
        "ci_test_results": ci_test_results,
        "z_test_results": z_test_results,
        "z_test_p_values": z_test_p_values,
        "mean_diff": mean_diff
    }
    if method:
        z_test_cor_res = check_dge_with_ztest(first_table, second_table, results=True, method=method)
        z_test_cor_p_val = check_dge_with_ztest(first_table, second_table, results=False, values=True, method=method)
        results["z_test_cor_res"] = z_test_cor_res
        results["z_test_cor_p_values"] = z_test_cor_p_val
        return pd.DataFrame(results)
    df_results = pd.DataFrame(results)
    df_results.to_csv(f"{save_results_table}.csv")


if __name__ == "__main__":
    while True:
        print("Enter the path to the file with the first cell type")
        first = str(input())
        print("Enter the path to the file with the second cell type")
        second = str(input())
        print("Enter the path to the results")
        res = str(input())
        print("Enter the method for multiple comparisons")
        method = str(input())
        create_table(first, second, res, method)
        print("You results are ready! Do you want to quit? y/n")
        command = str(input())
        command = command.lower()
        match command:
            case "yes" | "y":
                break
            case "no" | "n":
                pass
            case _:
                print("No command. Good bye")
                break
