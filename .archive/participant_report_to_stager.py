import argparse
import json
import requests
import logging
import sys
import pandas as pd


def get_parser():
    parser = argparse.ArgumentParser(
        "Takes a participant-wise report generated from 'participant_report_extraction.py', converts into a json and POSTs it to api/analyses/{analysis_id}/datasets/{dataset_id}/variants"
    )
    parser.add_argument(
        "--participant_report_path",
        type=str,
        required=True,
        help="path to participant-wise report csv",
    )
    parser.add_argument(
        "--analysis_id",
        type=int,
        required=True,
        help="Stager's analysis ID associated with the report",
    )
    parser.add_argument(
        "--dataset_participant_mapping",
        type=str,
        required=True,
        help="JSON mapping between participant codenames and dataset IDs, also from Stager",
    )
    parser.add_argument(
        "--user_id",
        type=int,
        required=True,
        help="User ID for Stager",
    )

    return parser


def replace_nan_with_none(record: dict[str, any]) -> dict[str, any]:
    record = {k: None if pd.isna(v) else v for k, v in record.items()}
    return record


def report_to_dict(csv_path: str) -> dict[str, any]:

    identifier_columns = ["analysis", "participant", "family"]

    variant_cols = [
        "chromosome",
        "position",
        "reference_allele",
        "alt_allele",
        "variation",
        "refseq_change",
        "depth",
        "conserved_in_20_mammals",
        "sift_score",
        "polyphen_score",
        "cadd_score",
        "gnomad_af",
        "ucsc_link",
        "gnomad_link",
        "gene",
        "info",
        "quality",
        "clinvar",
        "gnomad_af_popmax",
        "gnomad_ac",
        "gnomad_hom",
        "report_ensembl_gene_id",
        "ensembl_transcript_id",
        "aa_position",
        "exon",
        "protein_domains",
        "rsids",
        "gnomad_oe_lof_score",
        "gnomad_oe_mis_score",
        "exac_pli_score",
        "exac_prec_score",
        "exac_pnull_score",
        "spliceai_impact",
        "spliceai_score",
        "vest3_score",
        "revel_score",
        "gerp_score",
        "imprinting_status",
        "imprinting_expressed_allele",
        "pseudoautosomal",
        "number_of_callers",
        "old_multiallelic",
        "uce_100bp",
        "uce_200bp",
    ]

    variant_analysis_dataset_columns = [
        "zygosity",
        "burden",
        "alt_depths",
        "coverage",
        "genotype",
    ]

    ptp_df = pd.read_csv(csv_path)

    ptp_df["family"] = ptp_df["family"].apply(str)
    gt_dict = ptp_df[variant_analysis_dataset_columns].to_dict(orient="records")

    # replace nan/none's with null for consistency and mysql compatability
    for i in range(len(gt_dict)):
        gt_dict[i] = replace_nan_with_none(gt_dict[i])

    vt_dict = ptp_df[[col for col in variant_cols if col in ptp_df.columns]].to_dict(
        orient="records"
    )

    # replace nan/none's with null for consistency and mysql compatability
    for i in range(len(vt_dict)):
        vt_dict[i] = replace_nan_with_none(vt_dict[i])
        # nest genotype object inside its variant
        vt_dict[i]["genotype"] = gt_dict[i]

    identifier_dict = {k: ptp_df[k].iloc[0] for k in identifier_columns}

    all_dict = {**identifier_dict, "variant": vt_dict}

    return all_dict


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(asctime)s (line %(lineno)s): %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    args = get_parser().parse_args()

    logging.info(
        "Reading in {} and converting to json".format(args.participant_report_path)
    )
    report_dict = report_to_dict(args.participant_report_path)

    ptp_name = report_dict["participant"].split("_")[1]

    with open(args.dataset_participant_mapping) as json_f:
        mapping_dict = json.load(json_f)

    # placeholder
    dataset_id = mapping_dict.get(ptp_name)
    if not mapping_dict.get(ptp_name):
        logging.error("Participant codename not found in JSON mapping")
        sys.exit(1)

    logging.info(
        "Found Dataset ID '{}' for participant codename '{}'".format(
            dataset_id, mapping_dict.get(ptp_name)
        )
    )

    url = "http://localhost:5000/api/analyses/{}/datasets/{}/variants?user={}".format(
        args.analysis_id, dataset_id, args.user_id
    )
    logging.info("POSTing dictionary to endpoint: {}".format(url))
    response = requests.post(
        url=url,
        headers={"Content-Type": "application/json"},
        data=json.dumps(report_dict, indent=4),
    )

    try:
        response.raise_for_status()
        logging.info("Participant report successfully POSTed")
    except requests.exceptions.HTTPError as e:
        # get error msg from endpoint
        logging.error(
            "Error Code {} for {}: {} ".format(
                response.status_code, args.participant_report_path, response.json()
            )
        )
