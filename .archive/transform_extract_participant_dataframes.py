import os
import pandas as pd
import argparse
import sys
import numpy as np
from typing import Tuple, List, Iterator


def reshape_reports(
    samples: List[str], df: pd.DataFrame, report_fn: str
) -> Iterator[Tuple[str, str, str, pd.DataFrame]]:
    """
    creates participant-wise (wrt variants) dataframes from a cre report
    """

    for i, sample in enumerate(samples):

        sample_df = df.copy()

        # create a dictionary corresponding to the specific samples genotype columns
        genotype_cols = {
            gt.lower(): "{}.{}".format(gt, sample).lower()
            for gt in ["Zygosity", "Burden", "Alt_depths"]
        }

        # remove existing columns in the sample df corresponding to genotype columns
        sample_df.drop(
            list(df.filter(regex="zygosity|burden|alt_depths|gts|trio_coverage")),
            axis=1,
            inplace=True,
        )

        # iterate over each field, creating a placeholder in participant-wise dataframe, fetching the appropriate participant-genotype column from the original dataframe and assigning it to a stand alone column in the participant-wise dataframe
        for gt_field in genotype_cols:

            # create placeholder in sample df corresponding to zygosity, burden, and alt_depths per participant
            sample_df[gt_field] = np.nan

            # if field exists in columns, then assign to sample_df otherwise it will remain NA
            if genotype_cols[gt_field] in df.columns:
                # eg. sample_df[Zygosity] = df[zygosity.1389_ch0200]
                sample_df[gt_field] = df[genotype_cols[gt_field]]

        # retrieve appropriate GT location
        if "gts" in df.columns:
            gts = df["gts"].str.split(",", expand=True)

            if gts.shape[1] != len(samples):
                print(
                    "The number of extracted genotypes does not match the number of samples"
                )
                sys.exit(1)

            # assume the order of the genotypes is the same as the samples by indexing and replace 'gts'
            sample_df["genotype"] = gts[i].tolist()

        if "trio_coverage" in df.columns:
            trio_coverage = df["trio_coverage"].str.split("_", expand=True)

            if trio_coverage.shape[1] != len(samples):
                print(
                    "The number of extracted coverage fields does not match the number of samples"
                )
                sys.exit(1)

            # assume the order of the genotypes is the same as the samples by indexing and replace 'gts'
            sample_df["coverage"] = trio_coverage[i].tolist()

        # filter out homozygous reference variants or insufficent coverage variants from frequency
        sample_df = sample_df[
            (sample_df["zygosity"] != "-")
            & (sample_df["zygosity"] != "Insufficient coverage")
        ]

        # add identifying information based on the report
        sample_df["participant"] = sample
        family = os.path.basename(report_fn).split(".")[0]
        sample_df["family"] = family
        sample_df["analysis"] = os.path.basename(report_fn)

        # reorder columns without having to specify all columns in a list.
        # works in the case where the columns are not all found in the dataframe
        cols_to_move = [
            "participant",
            "family",
            "analysis",
            "chromosome",
            "position",
            "reference_allele",
            "alt_allele",
            "zygosity",
            "burden",
            "alt_depths",
            "genotype",
            "coverage",
        ]
        sample_df = sample_df[
            cols_to_move + [col for col in sample_df.columns if col not in cols_to_move]
        ]

        date_report_generated = os.path.basename(report_fn).split(".")[-2]

        yield sample, family, date_report_generated, sample_df


def preprocess_report(report: str) -> Tuple[List[str], pd.DataFrame]:
    """
    given a report path, reads in a as a dataframe and performs wrangling to normalize the dataframes
    this function was originally made to normalize the report dataframes such that they could be inserted into Stager's variant store database
    """

    # subsets report and returns samples using 'Zygosity.' columns

    if report.endswith(".csv"):
        sep = ","
    elif report.endswith(".tsv"):
        sep = "\t"
    else:
        print("Unknown fileformat and delimiter")
        sys.exit(1)

    try:
        df = pd.read_csv(report, sep=sep)
    except UnicodeDecodeError:
        print("UnicodeDecodeError on %s. Trying latin-1 decoding." % report)
        df = pd.read_csv(report, encoding="latin-1", sep=sep)
    except Exception as e:
        print(
            "Report '{}' could not be read in, please double check this is a valid csv!".format(
                report
            )
        )
        sys.exit(1)

    # print(df.shape)

    # these columns are 'wide' wrt the variants
    d = {"Zygosity": [], "Burden": [], "Alt_depths": []}

    # get all columns with Zygosity, Burden, or Alt_depths - these are unique for each participant
    for key in d:
        d[key].extend([col for col in df.columns if col.startswith(key)])

    # get sample names - preserved order for genotype and trio coverage
    samples = [col.replace("Zygosity.", "").strip() for col in d["Zygosity"]]

    # convert columns to lowercase
    df.columns = map(str.lower, df.columns)

    # rename ref and alt alleles to match model
    df = df.rename(
        columns={
            "ref": "reference_allele",
            "alt": "alt_allele",
            "ensembl_gene_id": "report_ensembl_gene_id",
        }
    )

    # convert variation values to lowercase
    df["variation"] = df["variation"].str.lower()

    # split 'position' into chromosome and position columns
    df[["chromosome", "position"]] = df["position"].str.split(":", expand=True)

    # make None's consistent (doesn't account for 0's though)
    df = df.replace({"None": None, np.nan: None})

    for col in ["conserved_in_20_mammals", "vest3_score", "revel_score", "gerp_score"]:
        if col in df.columns:
            df[col] = df[col].replace({".": None})

    df["depth"] = df["depth"].fillna(0)

    # weird column in that there are only Yes or NA's
    if "pseudoautosomal" in df.columns:
        df["pseudoautosomal"] = df["pseudoautosomal"].replace(
            {"Yes": 1, None: 0, np.nan: 0}
        )

    # remove hyperlink formatting
    for link_col in ["ucsc_link", "gnomad_link"]:

        # extract odd values b/w quotes
        try:
            df[link_col] = df[link_col].apply(lambda x: x.split('"')[1::2][0])
        # some columns don't actually have a link and therefore aren't splittable, eg. # results/14x/1473/1473.wes.2018-09-25.csv
        except IndexError as e:
            pass

    # take the last element after splitting on ';', won't affect non ';' delimited values
    df["clinvar"] = df["clinvar"].map(lambda x: str(x).split(";")[-1])

    # if coding persists. replace with appropriate value
    df["clinvar"] = df["clinvar"].replace(
        {
            "255": "other",
            "0": "uncertain",
            "1": "not-provided",
            "2": "benign",
            "3": "likely-benign",
            "4": "likely-pathogenic",
            "5": "pathogenic",
            "6": "drug-response",
            "7": "histocompatability",
        },
        regex=True,
    )
    # replace all forward slashes with '|' to ensure consistency
    df["clinvar"].replace({"\\/": "|"}, regex=True)

    # lower case
    df["clinvar"] = df["clinvar"].replace({"None": None, np.nan: None})
    df["clinvar"] = df["clinvar"].str.lower()

    # looks like report changed around 2021-01. before, spliceai_score contains | delimited impact and spliceai_impact contains a float
    # afterwards spliceai_score contains the float and spliceai_impact contains the | delimited score
    if "spliceai_score" in df.columns:
        # print(df["spliceai_score"])
        try:
            if any(df["spliceai_score"].str.contains("|")):
                df["spliceai_impact"] = df["spliceai_score"]
                df["spliceai_score"] = None
        # is a proper spliceai_score ie float. use a try catch since type checks don't work well (can be object, string or even sometimes numeric)
        # could force the dtype when reading in the csv but would prefer not to
        except AttributeError as e:
            pass

    return samples, df


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report_path",
        type=str,
        required=True,
        help="path to report csv",
    )

    return parser


if __name__ == "__main__":

    args = get_parser().parse_args()

    for arg in vars(args):
        print(arg + ": " + str(getattr(args, arg)))

    participants, df = preprocess_report(args.report_path)

    # which sample identifier do we want to use, right now 'sample' is actually {family}_{sample} resulting in the filename being
    # {family_sample}_{family}.csv
    # meaning 'family' is not actually required to be parsed out since it's a part of the 'sample' name
    for participant, family, date_report_generated, sample_df in reshape_reports(
        participants, df, args.report_path
    ):

        print("Processing sample {} from family {}".format(participant, family))
        # ptp_fn = '{}_{}.csv'.format(participant, family)
        ptp_fn = "{}_{}.csv".format(participant, date_report_generated)

        print("Saving to csv..")
        sample_df.to_csv(ptp_fn, index=False)
        print("Done saving csv.")

    print("Done")
