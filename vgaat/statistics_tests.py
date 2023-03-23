import pandas as pd
from mutations import CategorizedMutations
from variants import CategorizedVariantRecords
import utils
from collections import defaultdict


def mutation_site_count_vs_sample_type(
    cat_muts: CategorizedMutations,
    cat_vars: CategorizedVariantRecords,
    viral_strand,
    region,
    non_synonymous,
    allow_fishers,
):
    # Create list of earliest mutations
    earliest_mutations_category = "earliest mutations"
    cat_muts.pcintersect_category_groups_with_categories(
        earliest_mutations_category,
        "sample type",
        tuple(cat_muts.flattened_pcategory_group("patient earliest"))
    )

    working_category_name = earliest_mutations_category
    next_working_category_name = earliest_mutations_category
    if viral_strand != "all":
        next_working_category_name += f" - {viral_strand}"
        cat_muts.pcintersect_category_groups_with_categories(
            next_working_category_name,
            working_category_name,
            tuple(cat_muts.pget_category_group("viral lineage")[viral_strand]),
        )
        working_category_name = next_working_category_name
    if region != "all":
        next_working_category_name += f" - {region}"
        cat_muts.pcintersect_category_groups_with_categories(
            next_working_category_name,
            working_category_name,
            tuple(cat_muts.pget_category_group("regions")[region]),
        )
        working_category_name = next_working_category_name

    if non_synonymous:
        next_working_category_name += " - non-synonymous"
        cat_muts.pcintersect_category_groups_with_categories(
            next_working_category_name,
            working_category_name,
            tuple(cat_muts.pget_category("non-synonymous")),
        )
        working_category_name = next_working_category_name

    mutated_site_counts_per_sample_type = (
        cat_muts.mutation_sites_for_pccategory_group(working_category_name)
    )
    reference_site_counts_per_sample_type = dict(mutated_site_counts_per_sample_type)

    for key, value in reference_site_counts_per_sample_type.items():
        reference_site_counts_per_sample_type[key] = len(cat_muts.ref_sequence) - value

    conti_table_dict = {
        "Mutated": mutated_site_counts_per_sample_type,
        "Reference": reference_site_counts_per_sample_type,
    }

    contigency_test = pd.DataFrame(conti_table_dict)
    result = utils.row_pairwise_dependency_test(contigency_test, allow_fishers)

    return (
        "Is mutation site count dependent on sample type?",
        "Here, we will count the sites of mutation in the "
        "reference sequence and compare that number for each "
        "type of mutation. The following is the table we will "
        "we will use: ",
        utils.publish_results(result),
    )


# TODO Write individual nucleotide test
def individual_mutation_count_vs_sample_type(
    cat_muts: CategorizedMutations,
    cat_vars: CategorizedVariantRecords,
    viral_strand,
    region,
    non_synonymous,
    allow_fishers,
):
    results = []
    for location, location_based_mutations in cat_muts.pget_category_group("locations"):
        contingency_table_per_location = defaultdict(dict)
        for sample_type, sample_type_based_mutations in cat_muts.pget_category_group(
            "sample type"
        ):
            contingency_table_per_location[f"Mutations at {location}"][
                sample_type
            ] = len(
                cat_muts.intersect(
                    location_based_mutations, sample_type_based_mutations
                )
            )
            contingency_table_per_location["Reference"][sample_type] = len(
                location_based_mutations
            ) - len(cat_muts.pget_category("patients"))
        utils.publish_results(contingency_table_per_location, results)

    return "Location based mutations vs Sample type", results


# TODO Write distribution difference test (?)

tests = [
    (
        "Correlation between the number of mutation locations and sample types",
        mutation_site_count_vs_sample_type,
    )
]
