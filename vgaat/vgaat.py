#!/usr/bin/env python3

import os
import asyncio
import percache
import tempfile
import argparse
import logging
import variants
import mutations
import statistics_tests
import pickle
import shutil
import utils
import filters


async def main(args):
    logging.basicConfig(level=args.log.upper())
    logging.debug(f'Caching to "{args.cache_dir}".')
    if args.clear_cache:
        logging.info("Clearing previous cache...")
        shutil.rmtree(args.cache_dir)
        logging.info("Cache cleared")

    logging.info('Using call name "{0}"'.format(args.call_name))

    vcf_dir_path = os.path.abspath(args.vcf_dir)
    logging.info(f'Fetching VCF files from "{vcf_dir_path}"')
    lineage_file = os.path.abspath(args.lineage_file)
    logging.info(f'Fetching Lineage file from "{lineage_file}"')
    variant_organizer = variants.VariantRecordsOrganizer()
    variant_organizer.update(vcf_dir_path, lineage_file, args.sample_filename_re)
    logging.info("Building categorized variants...")
    categorized_variants = variant_organizer.build()
    logging.info("Done")
    mutation_organizer = mutations.MutationOrganizer(categorized_variants)

    logging.info(f"Using GenBank file from {args.ref_genbank}")
    mutation_organizer.update(args.ref_genbank, args.call_name)
    categorized_mutations_cache_path = os.path.join(
        args.cache_dir, os.path.basename(vcf_dir_path), "categorized_mutations.pickle"
    )
    if not os.path.exists(categorized_mutations_cache_path):
        logging.info("Building categorized mutations...")
        categorized_mutations = mutation_organizer.build()
        os.makedirs(os.path.dirname(categorized_mutations_cache_path))
        with open(categorized_mutations_cache_path, "wb") as fd:
            pickle.dump(categorized_mutations, fd)
    else:
        logging.info(
            f"Loading categorized mutations from {categorized_mutations_cache_path}"
        )
        with open(categorized_mutations_cache_path, "rb") as fd:
            categorized_mutations = pickle.load(fd)
    logging.info("Done")

    # TODO Add all categories as parameters
    # TODO How do we create a unanimous test suite???
    tester = utils.Tester(
        statistics_tests.tests,
        ["all", *categorized_mutations.pget_category_group("viral lineage").keys()],
        ["all", *categorized_mutations.pget_category_group("regions").keys()],
        [True, False],
        [not args.disable_fishers],
        categorized_mutations,
        categorized_variants,
        max_threads=args.threads,
    )
    results_cache_path = os.path.join(
        args.cache_dir, os.path.basename(vcf_dir_path), "results.pickle"
    )
    if not os.path.exists(results_cache_path):
        logging.info("Running all tests...")
        tester.run_all_async()
        results = tester.get_all_results()
        os.makedirs(args.cache_dir, exist_ok=True)
        with open(results_cache_path, "wb") as fd:
            pickle.dump(results, fd)
    else:
        logging.info(f"Loading test results from {results_cache_path}")
        with open(results_cache_path, "rb") as fd:
            results = pickle.load(fd)

    logging.info(f"Applying alpha filter of {args.alpha}")
    results = filters.filter_by_alpha(results, args.alpha)

    if not args.output:
        logging.debug("Outputting results to stdout...")
        utils.write_markdown_results(results)
    else:
        logging.debug(f'Outputting to "{args.output}"...')
        utils.write_markdown_results(results, md_results_path=args.output)
    logging.debug("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="VGAAT",
        description="Virus Genome Association Analytics Tools (VGAAT) \
            is a python program tool set containing a variety of associative \
                algorithms that may be run upon large amounts of VCFs.",
    )
    parser.add_argument(
        "vcf_dir", metavar="i", help="Path to directory containing VCF files"
    )
    parser.add_argument(
        "ref_genbank",
        metavar="a",
        help="The path to the NCBI GenBank file containing the reference used to \
            produce the VCF calls.",
    )
    parser.add_argument(
        "call_name",
        metavar="c",
        help="The call name to use when reading the VCF files.",
    )
    parser.add_argument(
        "lineage_file",
        metavar="l",
        help="The CSV file containing information on the samples lineage.",
    )
    parser.add_argument(
        "--sample-filename-re",
        metavar="-S",
        help="The regex used to interpret the individual sample filenames.",
        default=r"([BNE]{1,2})(\d+)(?:-D(\d+))?",
    )
    parser.add_argument(
        "--log",
        metavar="-L",
        help="Sets the verbosity of the program.",
        default="INFO",
    )
    parser.add_argument(
        "--cache-dir",
        metavar="-C",
        help="Set data cache location. Choose a persistent location if you'd like to \
            persist data after a run.",
        default="./tmp/VGAAT/data_cache",
    )
    parser.add_argument(
        "--disable-fishers",
        metavar="-X",
        help="Disables use of the Fisher's Exact Test even when it is possible.",
        default=False,
    )
    parser.add_argument(
        "--threads",
        metavar="-T",
        help="Number of threads to use when performing statistical tests.",
        default=16,
        type=int,
    )
    parser.add_argument(
        "--clear-cache", metavar="-S", help="Clears cache and then runs.", default=False
    )
    parser.add_argument(
        "--alpha",
        metavar="A",
        help="Filter results to be within given alpha value.",
        default=0.05,
    )
    parser.add_argument(
        "--output",
        metavar="-o",
        help="Where to output the results.",
    )
    # TODO Complete adding output and file format options

    args = parser.parse_args()
    try:
        cache = percache.Cache(os.path.join(tempfile.gettempdir(), "cache"))
    except PermissionError:
        cache = percache.Cache(args.cache_dir)
    asyncio.run(main(args))
