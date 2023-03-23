import hashlib
import os
import shutil
import sys
from scipy.stats import chi2_contingency, fisher_exact
import pandas as pd
import itertools
from multiprocessing.pool import ThreadPool
import logging


def accession_to_genbank_filename(accession):
    return accession.replace(".", "_") + ".gb"


def genbank_file_accession(gb_name):
    return gb_name.replace("_", ".").replace(".gb", "")


def row_pairwise_dependency_test(table: pd.DataFrame, allow_fishers=True):
    pairs = []
    indices = table.index
    results = {}
    for i in range(len(indices) - 1):
        for j in range(i + 1, len(indices)):
            pairs.append((indices[i], indices[j]))
    for a, b in pairs:
        row_pair_unmodified = table.loc[[a, b]]
        row_pair = row_pair_unmodified.loc[
            :, (row_pair_unmodified != 0).any(axis=0)
        ]  # Magic to drop columns that are empty...
        row_pair = row_pair.loc[(row_pair != 0).any(axis=1)]
        if row_pair.shape[0] == 2 and row_pair.shape[1] == 2 and allow_fishers:
            odds_ratio, p = fisher_exact(row_pair)
            results[f"{a} - {b}"] = {
                "type": "exact",
            }
        elif row_pair.shape[1] > 2 and row_pair.shape[0] >= 2:
            chi2, p, dof, expected = chi2_contingency(row_pair)
            results[f"{a} - {b}"] = {
                "type": "chi2",
            }
        elif row_pair.shape[1] == 0 and row_pair.shape[0] == 0:
            continue
        else:
            results[f"{a} - {b}"] = {
                "type": "chi2",
                "original_table": row_pair_unmodified,
            }
            p = "Error"
        results[f"{a} - {b}"]["p"] = p
        results[f"{a} - {b}"]["table"] = row_pair
    return results


# TODO Clean this up...
def publish_results(results: dict, results_list=None):
    results_list = results_list or []
    for result_info in results.values():
        results_list.append(result_info)
    return results_list


class Tester:
    def __init__(
        self,
        test_funcs,
        viral_strands,
        regions,
        synonymity,
        fishers,
        categorized_mutations,
        categorized_variants,
        max_threads=16,
    ):
        self.tests = test_funcs
        self.viral_strands = viral_strands
        self.regions = regions
        self.synonymity = synonymity
        self.fishers = fishers
        self.categorized_mutations = categorized_mutations
        self.categorized_variants = categorized_variants
        self._max_threads = max_threads
        self._results = {}

    def run_all_async(self):
        self._thread_pool = ThreadPool(processes=self._max_threads)
        param_product = itertools.product(
            self.tests, self.viral_strands, self.regions, self.synonymity, self.fishers
        )
        total = len(self.tests) * \
            len(self.viral_strands) * \
            len(self.regions) * \
            len(self.synonymity) * \
            len(self.fishers)
        runs = self._thread_pool.imap_unordered(
            self.run_test,
            param_product,
            chunksize=self._max_threads  # TODO Perhaps add more sophisticated param...
        )
        for running_test in runs:
            identifier, res = running_test
            self._results[identifier] = res
            logging.info(f"Test progress: {(len(self._results))/total:.1%}")

    def run_test(self, params):
        test_tuple, viral_strand, regions, synonymity, fishers = params
        test_alias, test_func = test_tuple
        logging.debug(
            "Running {0} with parameters: {1} {2} {3} {4}".format(
                test_alias, viral_strand, regions, synonymity, fishers
            )
        )
        res = test_func(
            self.categorized_mutations,
            self.categorized_variants,
            viral_strand,
            regions,
            synonymity,
            fishers
        )
        logging.debug("Completed running {0} with parameters: {1} {2} {3} {4}".format(
            test_alias, viral_strand, regions, synonymity, fishers
        ))
        return (test_alias, viral_strand, regions, synonymity, fishers), res

    def get_result_list(self, test_alias: str, viral_strand: str, regions,
                        synonymity: bool, fishers: bool):
        return self._results[test_alias, viral_strand, regions, synonymity, fishers]

    def get_all_results(self):
        return self._results


def write_markdown_results(result_groups: dict[str, list[dict]], md_results_path=None):
    if md_results_path:
        large_dir = md_results_path + ".large"
        if os.path.exists(large_dir):
            shutil.rmtree(large_dir)

    writer_p = sys.stdout
    if md_results_path:
        writer_p = open(md_results_path, "w")
    for test_id, result_list in result_groups.items():
        result_group_title, result_group_desc, results = result_list
        writer_p.write(f"# {result_group_title}\n\n")
        writer_p.write(f"{result_group_desc}\n\n")
        writer_p.write("Run with the following parameters:\n")
        writer_p.write(
            """ Test name: {0}; \
                Lineage: {1}; \
                Region: {2}; \
                Allow synonymous: {3}; \
                Allow Exact (Fishers): {4};\
                \n\n""".format(*test_id)
        )

        for result in results:
            writer_p.write(f"P-value: {result['p']:.5%}; \
                Test Type: {result['type']}\n\n")
            if result["table"].shape[0] < 10 and result["table"].shape[1] < 10:
                writer_p.write(f"{result['table'].to_markdown()}\n---\n\n")
            else:
                writer_p.write(
                    f"Table was too large {result['table'].shape} to display.\n"
                )
                if md_results_path:
                    large_table_dir = os.path.join(
                        md_results_path + ".large",
                        os.path.sep.join(map(str, test_id)).
                        replace(" ", "").
                        replace("\t", "").
                        replace("'", "").
                        replace("\"", "").
                        replace("(", "").
                        replace(")", "").
                        replace(",", "_")
                    )
                    filename = str(
                        hashlib.shake_128(
                            result["table"].to_markdown().encode("utf-8")
                        ).hexdigest(8)
                    )
                    filepath = os.path.join(
                        large_table_dir,
                        filename
                    )
                    same_num = 0
                    while os.path.exists(filepath + str(same_num) + ".csv"):
                        same_num += 1
                    large_table_path = filepath + str(same_num) + ".csv"
                    relative_table_path = large_table_path.replace(
                        os.path.abspath(md_results_path),
                        "." + os.path.sep + os.path.basename(md_results_path)
                    )
                    os.makedirs(large_table_dir, exist_ok=True)
                    result["table"].to_csv(path_or_buf=large_table_path)
                    writer_p.write(
                        "Table stored as CSV file. See at:\n"
                        f"[{relative_table_path}]({relative_table_path})"
                    )
                    writer_p.write("\n")
    if md_results_path:
        writer_p.close()
