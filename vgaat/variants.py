import csv
import re
import os
from typing import Iterable
import vcfpy
from collections import defaultdict
import logging


class VariantRecordsOrganizer:
    def __init__(self, on_change=None):
        self._on_change = on_change
        self._vcfs_dir = ""
        self._lineage_report_csv_path = ""
        self._vcf_filename_regex = ""

    def update(self, vcfs_dir, lineage_report_csv_path, vcf_filename_regex, **kwargs):
        self._vcfs_dir = vcfs_dir
        self._lineage_report_csv_path = lineage_report_csv_path
        self._vcf_filename_regex = vcf_filename_regex
        if self._on_change:
            self._on_change(self, **kwargs)

    def build(self):
        return CategorizedVariantRecords.from_files(
            self._vcfs_dir, self._lineage_report_csv_path, self._vcf_filename_regex
        )


class VariantRecordGroup:
    def __init__(self, patient_id, sample_type, sample_day, records: list):
        self._patient_id = patient_id
        self._sample_type = sample_type
        self._sample_day = sample_day
        self._variants = records

    @property
    def vcf_info(self):
        return self._patient_id, self._sample_type, self._sample_day

    @property
    def variants(self):
        return self._variants

    def __hash__(self):
        return hash(self.vcf_info)

    def __eq__(self, o):
        return self.vcf_info == o.vcf_info

    def __iter__(self):
        return iter(self.variants)


class CategorizedVariantRecords:
    def __init__(self, vcf_data, lineage_info):
        self._primitive_category_groups = {}
        self._primitive_category_groups["sample day"] = defaultdict(set)
        self._primitive_category_groups["patient id"] = defaultdict(set)
        self._primitive_category_groups["sample type"] = defaultdict(set)
        self._primitive_category_groups["viral lineage"] = defaultdict(set)
        self._all_vcf = vcf_data
        self._all_lineage_info = lineage_info
        for vcf_info, records in vcf_data.items():
            patient_id, sample_type, sample_day = vcf_info
            self._primitive_category_groups["sample day"][sample_day].add(records)
            self._primitive_category_groups["patient id"][patient_id].add(records)
            self._primitive_category_groups["sample type"][sample_type].add(records)
            if vcf_info in lineage_info:
                self._primitive_category_groups["viral lineage"][
                        lineage_info[vcf_info]
                    ].add(records)

    def get_category_groups(self, category: str):
        return self._primitive_category_groups[category]

    def get_all_vcf(self):
        return self._all_vcf

    def get_all_lineage_info(self):
        return self._all_lineage_info

    @staticmethod
    def intersect_each_with_one(each: dict, one: Iterable):
        result = {}
        for key, value in each.items():
            result[key] = value & one
        return result

    @staticmethod
    def get_sample_info_from_filename(
        vcf_filename, filename_regex=r"([BNE]{1,2})(\d+)(?:-D(\d+))?"
    ):
        filename_match = re.match(filename_regex, vcf_filename)
        if filename_match is None:
            logging.warn(
                f'The regex "{filename_regex}" did not match ' f'"{vcf_filename}"'
            )
            return None
        sample_type = filename_match.group(1)
        sample_day = (
            filename_match.group(3) or "1"
        )  # We assume it's day 1 in the case there were no annotations.
        sample_id = filename_match.group(2)
        return sample_id, sample_type, sample_day

    @staticmethod
    def read_lineage_csv(
        file_name: str, sample_info_column: str, viral_strand_column: str
    ):
        lineage_info = {}
        with open(file_name) as lineage_csv_fd:
            reader = csv.reader(lineage_csv_fd)
            header = None
            for row in reader:
                if header is None:
                    header = row
                    continue
                sample_info = CategorizedVariantRecords.get_sample_info_from_filename(
                    row[header.index("taxon")]
                )
                viral_strand_ver = row[header.index(viral_strand_column)]
                if viral_strand_ver:
                    lineage_info[sample_info] = viral_strand_ver
                else:
                    lineage_info[sample_info] = "unknown"
        return lineage_info

    @staticmethod
    def count_samples_in_dict_of_sets(dict_samples: dict[str, list]):
        result = {}
        for key, value in dict_samples.items():
            result[key] = len(value)
        return result

    @staticmethod
    def from_files(
        vcfs_path,
        lineage_report_csv_path,
        filename_regex=r"([BNE]{1,2})(\d+)(?:-D(\d+))?",
    ):
        lineage_dict = CategorizedVariantRecords.read_lineage_csv(
            lineage_report_csv_path, "taxon", "scorpio_call"
        )
        vcf_records = {}
        for vcf_filename in os.listdir(vcfs_path):
            vcf_full_path = os.path.abspath(os.path.join(vcfs_path, vcf_filename))
            if os.path.isdir(vcf_full_path):
                continue
            if os.path.splitext(vcf_full_path)[1] != ".vcf":
                continue
            sample_info = CategorizedVariantRecords.get_sample_info_from_filename(
                vcf_filename, filename_regex
            )
            if sample_info is None:
                continue

            patient_id, sample_type, sample_day = sample_info
            with open(vcf_full_path) as vcf_file_desc:
                curr_list = []
                for vcf_record in vcfpy.Reader(vcf_file_desc):
                    curr_list.append(vcf_record)
                vcf_records[sample_info] = VariantRecordGroup(*sample_info, curr_list)
        return CategorizedVariantRecords(vcf_records, lineage_dict)
