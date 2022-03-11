import argparse
import os
import json
import logging
import numpy as np
import urllib3

from datetime import datetime
from glob import glob
from PTQuery import *
from typing import Tuple, List, Iterator
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm


logging.basicConfig(
    filename=f"variant-upload-{datetime.today().strftime('%Y-%m-%d')}.log",
    filemode="w+",
    level=logging.INFO,
    format="[%(levelname)s] %(asctime)s (line %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def read_report_csv(report:str) -> pd.DataFrame:
    """
    read in a WES report and return a pandas dataframe
    """

    if report.endswith(".csv"):
        sep = ","
    elif report.endswith(".tsv"):
        sep = "\t"
    else:
        logging.error(f"Unknown fileformat and delimiter for {report}")
        raise ValueError(f"Unknown fileformat and delimiter for {report}")

    try:
        df = pd.read_csv(report, sep=sep)
    except UnicodeDecodeError:
        df = pd.read_csv(report, encoding="latin-1", sep=sep)
    except Exception as e: 
        logging.error(f"Could not read in {report}")
        raise ValueError(f"Report '{report}' could not be read in, please double check this is a valid csv!")

    return df 

def preprocess_report(report: str) -> Tuple[List[str], pd.DataFrame]:
    """
    given a report path, reads in as a dataframe and performs wrangling to normalize the dataframe
    this function was originally made to normalize the report dataframes such that they could be inserted into Stager's variant store database
    """

    df = read_report_csv(report)

    # these columns are 'wide' wrt the variants
    d = {"Zygosity": [], "Burden": [], "Alt_depths": []}

    # get all columns ending with Zygosity, Burden, or Alt_depths - these are wrt to each participant
    for key in d:
        d[key].extend([col for col in df.columns if col.startswith(key)])

    # get sample names - preserved order for genotype and trio coverage
    samples = [col.replace("Zygosity.", "").strip() for col in d["Zygosity"]]

    df.columns = map(str.lower, df.columns)

    # rename to match template
    df = df.rename(columns={"omim_gene_description": "omim_phenotype" })

    # convert variation values to lowercase
    df["variation"] = df["variation"].str.lower()

    # make None's consistent (doesn't account for 0's though)
    df = df.replace({"None": None, np.nan: None})

    # replace '0's 
    for bad_col in ['omim_phenotype', 'orphanet']:
        if bad_col in df.columns:
            df[bad_col] = df[bad_col].replace({'0' : None})

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

    # if coding persists. replace with appropriate value - credit to Madeline C.
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
    df["clinvar"] = df["clinvar"].replace({"\\/": "|"}, regex=True)
    df["clinvar"] = df["clinvar"].replace({"None": None, np.nan: None}).str.lower()

    # these should be ints/null
    if "number_of_callers" in df.columns:
        if "Number_of_callers" in df["number_of_callers"].astype(str).values:
            df["number_of_callers"] = None

    # these should be ints/null, see. 1615.wes.regular.2020-08-14.csv
    if "frequency_in_c4r" in df.columns:
        if "Frequency_in_C4R" in df["frequency_in_c4r"].astype(str).values:
            df["frequency_in_c4r"] = None

    # the former will eventually be deleted from the PT schema (meaning the report will fail if the col is present)
    for col in ["seen_in_c4r_samples", "c4r_wes_samples"]:
        if col in df.columns:
            df[col] = None

    # looks like report formatting changed around 2021-01. before, only spliceai_score is present and contains the impact delimited by '|'
    # afterwards, spliceai_score contains the float and spliceai_impact contains the | delimited score
    if ('spliceai_score' in df.columns) and ('spliceai_impact' not in df.columns):
        df["spliceai_impact"] = df["spliceai_score"]
        df["spliceai_score"] = None

    # check for duplicate variants wrt to pos, ref, and alt, these seem to affect a fraction of reports eg. 1666, 1743, and 1516
    dupes = df[df.duplicated(['position', 'ref', 'alt'], keep=False)].shape[0]

    if dupes:
        logging.info(f'Duplicate variants found for {report}')
        logging.info(f'# of variants before: {df.shape[0]}')
        df = df.drop_duplicates(['position', 'ref', 'alt'])
        logging.info(f'# of variants after: {df.shape[0]}')

    return samples, df


def reshape_reports(samples: List[str], df: pd.DataFrame, report_fn: str) -> Iterator[Tuple[str, str, str, pd.DataFrame]]:
    """
    creates participant-wise dataframes from a cre report
    may want to move this elsewhere since it could also be called for Stager report POSTing
    """

    for i, sample in enumerate(samples):

        sample_df = df.copy(deep = True)

        # create a dictionary corresponding to the specific samples genotype columns
        genotype_cols = {
            gt: "{}.{}".format(gt, sample).lower() for gt in ["zygosity", "burden", "alt_depths"]
        }

        # remove existing columns in the sample df corresponding to genotype columns
        sample_df.drop(list(df.filter(regex="zygosity|burden|alt_depths|trio_coverage")), axis=1, inplace=True)

        # iterate over each field, creating a placeholder in participant-wise dataframe, fetching the
        # appropriate participant-genotype column from the original dataframe and assigning it to a stand alone column in the participant-wise dataframe
        for gt_field, sample_specific_field in genotype_cols.items():

            # create placeholder in sample df corresponding to zygosity, burden, and alt_depths per participant
            sample_df[gt_field] = np.nan

            if genotype_cols[gt_field] in df.columns:
                # eg. sample_df[Zygosity] = df[zygosity.1389_ch0200]
                sample_df[gt_field] = df[sample_specific_field]

        # retrieve appropriate GT location
        if "gts" in df.columns:
            gts = df["gts"].str.split(",", expand=True)

            if gts.shape[1] != len(samples):
                raise IndexError("The number of extracted genotypes does not match the number of samples")

            # assume the order of the genotypes is the same as the samples by indexing and replace 'gts' ;misnomer
            sample_df["gts"] = gts[i]


        if "trio_coverage" in df.columns:
            trio_coverage = df["trio_coverage"].astype(str)

            # replace dates, eg. 2003-05-06, 2004-10-03..... these are converted to dates by excel from eg. 03-05-06, 04-10-03
            # checks for values in this column that also fit eg. 11/34/54, that way we know the former was improperly converted to dates
            if any(trio_coverage.str.contains('[0-9]{4}-[0-9]{2}-[0-9]{2}', regex= True, na=False)) and any(trio_coverage.str.contains("\\/")):
        
                trio_coverage = trio_coverage.apply(lambda x: x[2:] if x.startswith('20') and '-' in x else x)
                trio_coverage = trio_coverage.str.replace('-', '/')
        
            if any(trio_coverage.str.contains("_")):
                trio_coverage = trio_coverage.str.split("_", expand=True)

                if trio_coverage.shape[1] != len(samples):
                    raise IndexError("The number of extracted coverage fields does not match the number of samples")

                # assume the order of the genotypes is the same as the samples by indexing and replace 'trio_coverage', misnomer
                sample_df["trio_coverage"] = trio_coverage[i].astype(int)

            elif any(trio_coverage.str.contains("\\/")):
                trio_coverage = trio_coverage.str.split("\\/", expand=True)

                if trio_coverage.shape[1] != len(samples):
                    raise IndexError("The number of extracted coverage fields does not match the number of samples")

                sample_df["trio_coverage"] = trio_coverage[i].astype(int)

            else:
                # singleton
                sample_df["trio_coverage"] = trio_coverage

        # filter out homozygous reference variants or insufficent coverage variants from frequency
        sample_df = sample_df[
            (sample_df["zygosity"] != "-") & (sample_df["zygosity"] != "Insufficient coverage") & (sample_df["zygosity"] != "Insufficient_coverage")
        ]

        # replace x/y/mt invalid zygosity values
        sample_df.loc[ ~sample_df["zygosity"].isin(["Het", "Hom", None]), "zygosity" ] = None

        # remove MT variants
        sample_df = sample_df[~sample_df["position"].astype(str).str.startswith("MT")]

        # add identifying information based on the report
        family = os.path.basename(report_fn).split(".")[0]

        # 'gts' and 'trio_coverage' should be 'gt' and 'coverage', respectively but are retained to be in line with the example PT singleton report 
        cols_to_move = ["position", "ref", "alt", "zygosity", "burden", "alt_depths", "gts", "trio_coverage"]

        sample_df = sample_df[cols_to_move + [col for col in sample_df.columns if col not in cols_to_move]]

        date_report_generated = os.path.basename(report_fn).split(".")[-2]

        yield sample, family, date_report_generated, sample_df


def format_for_phenotips(report: pd.DataFrame) -> pd.DataFrame:
    """
    Formatting for phenotips, specifically presence of columns and adding a header
    Credits to Conor Klamann for original code.
    """

    report.rename({"position": "#position"}, axis=1, inplace=True)

    missing_cols = set(ALLOWED_FIELDS).difference(set(report.columns.values))

    for col in missing_cols:
        report[col] = None

    extra_cols = set(report.columns.values).difference(set(ALLOWED_FIELDS))

    report.drop(columns=extra_cols, inplace=True)

    return extra_cols, missing_cols, report


def json_to_df_logs(master_list: list) -> pd.DataFrame:
    """
    convert list of dictionaries/json to a tidy Dataframe for analytics
    """
    df = pd.json_normalize(master_list)
    df = df.explode("participants")
    merged_df = pd.concat([df.drop("participants", axis=1), df["participants"].apply(pd.Series)], axis=1)
    return merged_df


def create_conn(authentication_method,
                username: Optional[str],
                password: Optional[str],
                data=Optional[dict]) -> PTQuery:
    """
    creates and returns an instance of PTQuery depending on authentication credentials
    """

    # staging instance
    if authentication_method == "basic":
        if not args.username and not args.password:
            raise ValueError( "Please make sure both username and password is passed in when using basic authentication")
        else:
            BASE_REQUEST_ARGS = {
                "headers": {},
                "verify": False,
            }

            base_url = "https://dev.phenotips.genomics4rd.ca"

            query = PTQuery(base_url=base_url, base_request_args=BASE_REQUEST_ARGS, username=username, password=password)

    elif authentication_method == "auth0":

        BASE_REQUEST_ARGS = {}
        base_url = "https://phenotips.genomics4rd.ca"
        bearer_token = get_bearer_token(data)

        query = PTQuery(base_url=base_url, base_request_args=BASE_REQUEST_ARGS, bearer_token=bearer_token)

    return query

def get_patient_id_from_mapping(family_participant_id:str, df, col_to_check = 'external_id') -> str:
    """
    query the merged df (intersection b/w family_reports and 
    Magda's curated phenotips list) with sample name to get internal Phenotips identifier
    """
    if family_participant_id in df[col_to_check].tolist():
        pt_id = df[df[col_to_check] == family_participant_id]['report_id'].values[0]
        return pt_id
    else:
        return None

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--username", type=str, required=False, help="basic auth username" )
    parser.add_argument("--password", type=str, required=False, help="basic auth password")
    parser.add_argument("--report-dir", type=str, required=False, help="path to folder of reports" )
    parser.add_argument("--report-path", type=str, required=False, help="path to a single report")
    parser.add_argument("--authentication-method", type=str, required=True, help="Basic or Auth0")
    parser.add_argument("--auth0-token-data", required=False, help="body required for auth0 authentication" )
    parser.add_argument("--resume-from-json", type=str, required=False, help="skips external ids found in this json")
    parser.add_argument("--mapping-file", type=str, required=False, help="use mapping dataframe instead of querying Phenotips for internal identifiers")
    return parser

if __name__ == "__main__":

    args = get_parser().parse_args()

    auth0_data = (
        json.loads(args.auth0_token_data) if args.auth0_token_data is not None else None
    )

    query = create_conn(
        authentication_method=args.authentication_method, 
        username=args.username, 
        password=args.password,
        data=auth0_data,
    )

    logging.info(f"Successfully authenticated as {args.username if args.username is not None else auth0_data['username']}")

    if not args.report_dir and not args.report_path:
        logging.error("Either 'report-dir' or 'report-path' arguments must be passed in.")
        raise ValueError("Either 'report-dir' or 'report-path' arguments must be passed in.")
    if args.report_dir:
        logging.info(f"Report directory passed in: {args.report_dir}")
        fn = f"variant-store-results-{datetime.today().strftime('%Y-%m-%d')}"
        report_dir = glob(os.path.join(args.report_dir, "**/*csv"), recursive=True)
    elif args.report_path:
        logging.info(f"Report path passed in: {args.report_path}")
        fn = f"variant-store-results-{os.path.basename(os.path.splitext(args.report_path)[0])}-on-{datetime.today().strftime('%Y-%m-%d')}"
        report_dir = [args.report_path]

    master_list = []

    folder_to_store_csvs = "demultiplexed_reports"

    if not os.path.exists(folder_to_store_csvs):
        os.makedirs(folder_to_store_csvs, exist_ok=True)

    if args.resume_from_json:
        logging.info(f'Resuming report pre-processing and POSTing from {args.resume_from_json}')
        with open(args.resume_from_json) as f:
            master_list = json.load(f)
        resume_df = json_to_df_logs(master_list)

    if args.mapping_file:
        logging.info(f'Using mapping file {args.mapping_file} to obtain phenotip identifiers...')
        mapping_df = pd.read_csv(args.mapping_file)

    with logging_redirect_tqdm():
        for report in tqdm(report_dir):

            if args.resume_from_json:
                if report in resume_df["report_name"].tolist():
                    logging.info(f"Skipping {report}...")
                    continue

            report_dict = {"report_name": None, "family": None, "participants": []}
            report_dict["report_name"] = report

            participants, df = preprocess_report(report)

            logging.info(f"Report Name: {report}")

            # create participant-wise df from the normalized report and POST
            for (family_participant_identifier, family,  date_report_generated, sample_df) in reshape_reports(participants, df, report):

                ptp_dict = {}

                report_dict["family"] = family

                if not args.mapping_file:
                    pt_id = str(query.get_internal_id_by_external_id(family_participant_identifier))
                else:
                    pt_id = get_patient_id_from_mapping(family_participant_identifier, mapping_df, col_to_check = 'external_id')

                ptp_dict["eid"] = family_participant_identifier
                ptp_dict["iid"] = pt_id

                if not pt_id:
                    logging.info(f"No Phenotips identifier found for {family_participant_identifier}")
                    ptp_dict["variants_found"] = None
                    ptp_dict["missing_cols"] = None
                    ptp_dict["extra_cols"] = None
                    ptp_dict["post_status_code"] = None

                elif pt_id.startswith("P"):

                    variants_exist = 0 #query.get_variant_info(pt_id)

                    if variants_exist:
                        logging.info(f"Variants found for {family_participant_identifier}")
                        ptp_dict["variants_found"] = 1

                    # no variants found so report will be formatted and POSTed
                    else:

                        ptp_dict["variants_found"] = 0

                        extra_cols, missing_cols, formatted_ptp_report = format_for_phenotips(sample_df)

                        ptp_dict["missing_cols"] = list(missing_cols)
                        ptp_dict["extra_cols"] = list(extra_cols)

                        ptp_fn = os.path.join(
                            folder_to_store_csvs,
                            f"{family_participant_identifier}_{date_report_generated}-formatted.csv",
                        )

                        formatted_ptp_report.to_csv(ptp_fn, index=False)

                        post_status_code = query.clean_and_post_report(pt_id, ptp_fn)

                        if post_status_code != 200:
                            logging.error(f"Report POST failed for {family_participant_identifier} with code {post_status_code}")
                        ptp_dict["post_status_code"] = post_status_code

                report_dict["participants"].append(ptp_dict)

            master_list.append(report_dict)

            # re-dump json with each report, is there a better way to do this? eg. append
            with open(f"{fn}.json", "w") as f:
                json.dump(master_list, f, indent=4)

        json_to_df_logs(master_list).to_csv(f"{fn}.csv", index=False)

    logging.info("Done!")

