from collections import defaultdict
from typing import Callable, Union
from Bio.Seq import MutableSeq, translate
from Bio import GenBank
from frozendict import frozendict
from utils import accession_to_genbank_filename
import logging
from functools import cache


class Mutation:
    def __init__(
        self,
        vcf_record,
        sample_id,
        sample_type,
        sample_day,
        viral_lineage,
        alt,
        sample_call_name="unknown",
    ):
        self._sample_viral_lineage = viral_lineage
        self._patient_id = sample_id
        self._sample_type = sample_type
        self._sample_day = sample_day
        self._start = vcf_record.begin
        self._end = vcf_record.end or self.start
        self._reference = vcf_record.REF
        self._call = vcf_record.call_for_sample[sample_call_name]
        self._genotype = self.call.gt_bases[0]  # Haploid
        self._all_alts = vcf_record.ALT
        self._alt = alt.value
        self._ref_obs = self.call.data["RO"]
        self._depth = self.call.data["DP"]

        if len(self.genotype) > len(self.ref):
            self._variant_type = "insertion"
        elif len(self.genotype) < len(self.ref):
            self._variant_type = "deletion"
        elif len(self.genotype) == len(self.ref):
            if len(self.genotype) == 1:
                self._variant_type = "SNV"
            else:
                self._variant_type = "MNV"

        self._hash = hash(
            (
                self._patient_id,
                self._sample_type,
                self._sample_day,
                self._start,
                self._end,
                self._reference,
                self._alt,
            )
        )

    @property
    def patient_id(self):
        return self._patient_id

    @property
    def sample_type(self):
        return self._sample_type

    @property
    def sample_day(self):
        return self._sample_day

    @property
    def sample_viral_lineage(self):
        return self._sample_viral_lineage

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def ref(self):
        return self._reference

    @property
    def call(self):
        return self._call

    @property
    def genotype(self):
        return self._genotype

    @property
    def alternatives(self):
        return self._all_alts

    @property
    def ref_obs(self):
        return self._ref_obs

    @property
    def depth(self):
        return self._depth

    @property
    def variant_type(self):
        return self._variant_type

    def __len__(self):
        return self.depth

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return (
            self.patient_id == other.get_patient_id
            and self.sample_type == other.get_sample_type
            and self.sample_day == other.sample_day
            and self.start == other.start
            and self.end == other.end
            and self.ref == other.reference()
            and self.call == other.call
            and self.genotype == other.genotype
            and self.alternatives == other.alternatives
            and self.ref_obs == other.ref_obs
            and self.depth == other.depth
            and self.variant_type == other.variant_type
            and self._alt == other._alt_nucleotide
        )

    @staticmethod
    def vcf_reads_to_mutations(
        vcf_reads: dict, accession, sample_lineage, sample_call_name="unknown"
    ):
        for sample_info, vcf_records in vcf_reads.items():
            sample_id, sample_type, sample_day = sample_info
            for vcf_record in vcf_records:
                if vcf_record.CHROM != accession:
                    continue
                for alternate in vcf_record.ALT:
                    yield Mutation(
                        vcf_record,
                        sample_id,
                        sample_type,
                        sample_day,
                        sample_lineage[sample_info]
                        if sample_info in sample_lineage
                        else None,
                        alternate,
                        sample_call_name,
                    )


class CategorizedMutations:
    def __init__(self, genbank_record, mutations: list[Mutation]):
        self._genbank_record = genbank_record
        self._primitive_category_groups = {}
        self._primitive_categories = {}
        self._custom_category_groups = {}
        self._custom_categories = {}
        genes_and_utrs_locs = {}

        logging.debug("Categorized mutations being constructed.")
        logging.debug("Cataloguing features...")
        for genbank_feature in self.get_genbank_record().features:
            f_key = genbank_feature.key
            if f_key == "gene" or f_key.endswith("UTR"):
                f_qual = f_key
                f_start, f_end = genbank_feature.location.split("..")
                f_start = int(f_start)
                f_end = int(f_end)
                for qualifier in genbank_feature.qualifiers:
                    if qualifier.key == "/gene=":
                        f_qual = qualifier.value
                        break
                if (f_start, f_end) in genes_and_utrs_locs:
                    raise Exception("Feature with duplicate location detected.")
                genes_and_utrs_locs[(f_start, f_end)] = f_qual
        logging.debug("Done")

        self._accession = genbank_record.accession

        self._primitive_category_groups["locations"] = defaultdict(set)
        self._primitive_category_groups["allele change"] = defaultdict(set)
        self._primitive_category_groups["sample type"] = defaultdict(set)
        self._primitive_category_groups["sample day"] = defaultdict(set)
        self._primitive_category_groups["variant type"] = defaultdict(set)
        self._primitive_category_groups["patient"] = defaultdict(set)
        self._primitive_category_groups["viral lineage"] = defaultdict(set)
        self._primitive_category_groups["regions"] = defaultdict(set)
        self._primitive_category_groups["patient earliest"] = defaultdict(set)

        self._primitive_categories["non-synonymous"] = set()
        self._primitive_categories["all"] = set(mutations)

        logging.debug("Organizing mutations...")
        num_of_mutations = len(mutations)
        done = 0
        for mutation in mutations:
            self._primitive_category_groups["locations"][
                (mutation.start, mutation.end)
            ].add(mutation)
            self._primitive_category_groups["allele change"][
                (mutation.ref, mutation.genotype)
            ].add(mutation)
            self._primitive_category_groups["sample type"][mutation.sample_type].add(
                mutation
            )
            self._primitive_category_groups["variant type"][mutation.variant_type].add(
                mutation
            )
            self._primitive_category_groups["sample day"][mutation.sample_day].add(
                mutation
            )
            if mutation.sample_viral_lineage:
                self._primitive_category_groups["viral lineage"][
                    mutation.sample_viral_lineage
                ].add(mutation)
            # Maintain all samples having the same day and by the end, that day
            # is the latest... If we start with nothing, then everything in the
            # set has the same day. This is base case. Then, when we add
            # something to the set, if one of the previous set has a older time
            # stamp, they all do (I.H). We clear the set and add the more
            # recent version. If the new timestamp is older then any arbitrary
            # item in the current set, then it's older than all items in the
            # set by I.H. We do not add this item. The last case if the item's
            # timestamp is the same as any arbitrary one. Therefore, all
            # items in the set will be the same. Since we are replacing
            # all older items with newer ones, the last update should contain
            # the latest set.
            if self._primitive_category_groups["patient earliest"][mutation.patient_id]:
                earliest_day = next(
                    iter(
                        self._primitive_category_groups["patient earliest"][
                            mutation.patient_id
                        ]
                    )
                ).sample_day
            else:
                earliest_day = None
            if earliest_day is None or earliest_day > mutation.sample_day:
                self._primitive_category_groups["patient earliest"][
                    mutation.patient_id
                ].clear()
            if earliest_day is None or mutation.sample_day == earliest_day:
                self._primitive_category_groups["patient earliest"][
                    mutation.patient_id
                ].add(mutation)
            self._primitive_category_groups["patient"][mutation.patient_id].add(
                mutation
            )
            mutable_sequence = MutableSeq(self.ref_sequence)
            for relative_mutation, mutation_position in enumerate(
                range(mutation.start, mutation.end + 1)
            ):
                mutable_sequence[mutation_position] = mutation.genotype[
                    relative_mutation
                ]
            if translate(mutable_sequence) != translate(self.ref_sequence):
                self._primitive_categories["non-synonymous"].add(mutation)

            for location in genes_and_utrs_locs.keys():
                start, end = location
                if start <= mutation.start and end >= mutation.end:
                    self._primitive_category_groups["regions"][
                        (genes_and_utrs_locs[location], location)
                    ].add(mutation)
            done += 1
            logging.info(
                f"Organized {done}/{num_of_mutations} \
                    ({done/num_of_mutations:.2%}) mutations"
            )

        self._primitive_categories["non-synonymous"] = set()
        self._primitive_categories["all"] = set(mutations)

        self._primitive_category_groups = recursive_freeze_dict(
            self._primitive_category_groups
        )

    def get_filename(self):
        """
        Returns the accession filename.
        """
        return accession_to_genbank_filename(self._accession)

    def get_accession(self):
        return self._accession

    def get_genbank_record(self):
        return self._genbank_record

    @property
    def ref_sequence(self):
        return self._genbank_record.sequence

    def pget_category_or_category_group(self, property_name: str):
        return (
            self._primitive_category_groups[property_name]
            or self._primitive_categories[property_name]
        )

    def pget_category(self, property_name: str):
        return self._primitive_categories[property_name]

    def pget_category_group(self, property_name: str):
        return self._primitive_category_groups[property_name]

    @property
    def all(self):
        return self._primitive_category_groups["all"]

    @property
    def total_mutations(self):
        return len(self.all())

    @cache
    def mutations_per_site_for_pcategory_group(self, category_group: str):
        categorized_mutations = self.pget_category_group(category_group)
        if isinstance(categorized_mutations, dict):
            resulting_dict = {}
            for key, value in categorized_mutations.items():
                resulting_dict[key] = self.count_mutations_per_site(value)
            return resulting_dict
        else:
            return self.count_mutations_per_site(categorized_mutations)

    @cache
    def flattened_pcategory_group(self, category: str) -> set:
        categorized_mutations = self.pget_category_group(category)
        flattened = set()
        for key, set_value in categorized_mutations.items():
            flattened.update(set_value)
        return flattened

    @cache
    def mutation_sites_for_pccategory_group(self, category: str):
        categorized_mutations = self.pcget_or_create_category_groups(category)
        resulting_dict = {}
        for key, value in categorized_mutations.items():
            resulting_dict[key] = CategorizedMutations.count_mutation_sites(value)
        return resulting_dict

    @cache
    def mutations_for_pccategory_group(
        self, category: str, assembler: Union[Callable, None] = None, **kwargs: dict
    ):
        dictionary_set = self.pcget_or_create_category_groups(category)
        return CategorizedMutations.count_mutations_for_category_group(dictionary_set)

    @cache
    def punion_category_groups_with_category(self, category_groups: str, category: str):
        category_set = set(self.pget_category(category))
        each = self.pget_category_group(category_groups)
        one = category_set
        result = {}
        for subcat, subcat_mutations in each.items():
            result[subcat] = subcat_mutations.union(one)
        return result

    @cache
    def pintersect_category_groups_with_category(
        self, category_groups: str, category: str
    ):
        category_set = set(self.pget_category(category))
        each = self.pget_category_group(category_groups)
        result = {}
        for subcat, subcat_mutations in each.items():
            result[subcat] = subcat_mutations.intersection(category_set)
        return result

    def pcget_or_create_category_groups(
        self, category_group: str, assembler: Union[Callable, None] = None, **kwargs
    ):
        if category_group in self._primitive_category_groups:
            return self.pget_category_group(category_group)
        if category_group not in self._custom_category_groups:
            if not assembler:
                raise Exception(
                    "Assembler was not provided for a non-existent custom category "
                    f'groups "{category_group}".'
                )
            result = assembler(**kwargs)
            if not isinstance(result, dict):
                raise Exception("Assembler result was not a group of categories...")
            self._custom_category_groups[category_group] = result
        return self._custom_category_groups[category_group]

    @cache
    def pcunion_category_groups_with_categories(
        self,
        union_name: str,
        category_group: str,
        category: tuple[Mutation],
        assembler: Union[Callable, None] = None,
        assembler_args: dict = {},
    ):
        category_set = set(category)
        if union_name in self._custom_category_groups:
            return self._custom_category_groups[union_name]

        subcategories = self.pcget_or_create_category_groups(
            category_group, assembler, **assembler_args
        )
        result = {}
        for subcat, mutations in subcategories.items():
            result[subcat] = mutations.union(category_set)

        self._custom_category_groups[union_name] = result
        return result

    @cache
    def pcintersect_category_groups_with_categories(
        self,
        intersect_name: str,
        category_group: str,
        category: tuple[Mutation],
        assembler: Union[Callable, None] = None,
        assembler_args: dict = {},
    ):
        category_set = set(category)

        if intersect_name in self._custom_category_groups:
            return self._custom_category_groups[intersect_name]

        subcategories = self.pcget_or_create_category_groups(
            category_group, assembler, **assembler_args
        )
        result = {}
        for subcat, mutations in subcategories.items():
            result[subcat] = mutations.intersection(category_set)
        self._custom_category_groups[intersect_name] = result
        return result

    @cache
    def pcbuild_contingency(
        self, category_group_a_name: str, category_group_b_name: str, modifier: Callable
    ):
        category_group_a = self.pcget_or_create_category_groups(category_group_a_name)
        category_group_b = self.pcget_or_create_category_groups(category_group_b_name)
        proto_contingency = defaultdict(dict)
        for a_group, mutations_cat_group_a in category_group_a.items():
            for b_group, mutations_cat_group_b in category_group_b.items():
                proto_contingency[a_group][b_group] = modifier(
                    self.intersect(
                        tuple(mutations_cat_group_a), tuple(mutations_cat_group_b)
                    )
                )
        return proto_contingency

    def clear_custom_categories(self):
        self._custom_category_groups.clear()

    @staticmethod
    @cache
    def count_mutations_per_site(mutations: tuple[Mutation]):
        mutation_sites = defaultdict(int)
        for mutation in mutations:
            for position in range(mutation.start, mutation.end + 1):
                mutation_sites[position] += 1
        return mutation_sites

    @staticmethod
    @cache
    def count_mutation_sites(mutations: tuple[Mutation]):
        mutation_sites = set()
        for mutation in mutations:
            for nucleotide_position in range(mutation.start, mutation.end + 1):
                mutation_sites.add(nucleotide_position)
        return len(mutation_sites)

    @staticmethod
    @cache
    def intersect(mutation: tuple[Mutation], *mutations: tuple[tuple[Mutation]]):
        intersect = set(mutation)
        for add_mut in mutations:
            intersect = intersect & set(add_mut)
        return intersect

    @staticmethod
    @cache
    def union(mutation: tuple[Mutation], *mutations: tuple[tuple[Mutation]]):
        union = set(mutation)
        for add_mut in mutations:
            union = union | set(add_mut)
        return union

    @staticmethod
    def count_mutations_for_category_group(
        group_of_mutations: dict[str, set[Mutation]]
    ):
        resulting_dict = {}
        for key, value in group_of_mutations.items():
            resulting_dict[key] = len(value)
        return resulting_dict

    @staticmethod
    def build_with_args(
        categorized_vcf_records, genbank_file, sample_call_name="sample_call_name"
    ):
        genbank_record = None
        with open(genbank_file) as genbank_fd:
            genbank_record = next(
                GenBank.parse(genbank_fd)
            )  # Assume there's only one accession per file.

        return CategorizedMutations(
            genbank_record,
            list(
                Mutation.vcf_reads_to_mutations(
                    categorized_vcf_records.get_all_vcf(),
                    genbank_record.version,  # type: ignore
                    categorized_vcf_records.get_all_lineage_info(),
                    sample_call_name,
                )
            ),
        )


class MutationOrganizer:
    def __init__(self, categorized_vcf_records, on_change=None):
        self.genbank_path = "unknown"
        self.sample_call_name = "unknown"
        self._on_change = on_change
        self._categorized_vcf_records = categorized_vcf_records

    def update(self, genbank_path, sample, **kwargs):
        self.genbank_path = genbank_path
        self.sample_call_name = sample
        if self._on_change:
            self._on_change(self, **kwargs)

    def build(self):
        return CategorizedMutations.build_with_args(
            self._categorized_vcf_records, self.genbank_path, self.sample_call_name
        )


def recursive_freeze_dict(to_freeze: dict) -> frozendict:
    for key, value in to_freeze.items():
        if isinstance(value, dict):
            to_freeze[key] = recursive_freeze_dict(value)
        elif isinstance(value, list):
            to_freeze[key] = tuple(value)
    return frozendict(to_freeze)


def mutations_collection_to_set(
    data: Union[tuple[Mutation], set[Mutation]]
) -> set[Mutation]:
    return set(data) if isinstance(data, tuple) else data
