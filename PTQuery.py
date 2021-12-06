from json import dumps
from os import path
from re import sub

import requests
import pandas as pd

""" headers, in order, from template singleton report """
ALLOWED_FIELDS = [
    "#Position",
    "UCSC_Link",
    "GNOMAD_Link",
    "Ref",
    "Alt",
    "Zygosity",
    "Gene",
    "Burden",
    "gts",
    "Variation",
    "Info",
    "Refseq_change",
    "Depth",
    "Quality",
    "Alt_depths",
    "Trio_coverage",
    "Ensembl_gene_id",
    "Gene_description",
    "omim_phenotype",
    "omim_inheritance",
    "Orphanet",
    "Clinvar",
    "Frequency_in_C4R",
    "Seen_in_C4R_samples",
    "HGMD_id",
    "HGMD_gene",
    "HGMD_tag",
    "HGMD_ref",
    "Gnomad_af_popmax",
    "Gnomad_af",
    "Gnomad_ac",
    "Gnomad_hom",
    "Ensembl_transcript_id",
    "AA_position",
    "Exon",
    "Protein_domains",
    "rsIDs",
    "Gnomad_oe_lof_score",
    "Gnomad_oe_mis_score",
    "Exac_pli_score",
    "Exac_prec_score",
    "Exac_pnull_score",
    "Conserved_in_20_mammals",
    "SpliceAI_impact",
    "SpliceAI_score",
    "Sift_score",
    "Polyphen_score",
    "Cadd_score",
    "Vest3_score",
    "Revel_score",
    "Gerp_score",
    "Imprinting_status",
    "Imprinting_expressed_allele",
    "Pseudoautosomal",
    "Number_of_callers",
    "Old_multiallelic",
    "UCE_100bp",
    "UCE_200bp",
]


class PTQuery:
    """
    Simple class for making requests to the PhenoTips API

    Attributes:
        base_url: The url of the PT endpoint, should end with top-level domain, no slash. E.g., phenotips.example.ca
        base_request_args: dict of kwargs to pass to the requests library. Should at minimum include {"headers": {"Authorization": <...>}}.
    """

    def __init__(self, base_url: str, base_request_args: dict):
        self.base_url = base_url
        self.base_request_args = base_request_args

    def clean_and_post_report(self, patient_id: str, report_path: str) -> str:
        """clean and post report for a patient, prints response status code"""
        cleaned_path = clean_report(report_path)
        filename = path.basename(cleaned_path)
        args = {
            "data": {
                "metadata": dumps(
                    {
                        "patientId": patient_id,
                        "refGenome": "GRCh37",
                        "fileName": filename,
                    }
                ),
            },
            "files": {"fileStream": (None, open(cleaned_path, "rb"))},
            "timeout": 30,
            **self.base_request_args,
        }
        res = requests.put(
            f"{self.base_url}/rest/variant-source-files/patients/{patient_id}/files/{filename}",
            **args,
        )
        return res.status_code

    def delete_report(self, patient_id: str, filename: str) -> str:
        """remove a report from a patient record (usually b/c of error when uploading)"""
        res = requests.delete(
            f"{self.base_url}/rest/variant-source-files/patients/{patient_id}/files/{filename}",
            **self.base_request_args,
        )
        return res.status_code

    def get_latest_job_metadata_for_patient(
        self, patient_id: str, complete=False, params={}
    ) -> str:
        """search metadata for participant's latest upload attempt"""

        PAGE_SIZE = 25

        if not params:
            params = {"patientLimit": PAGE_SIZE, "patientOffset": 0}
        if complete:
            params["procStatus"] = "COMPLETE"

        res = requests.get(
            "http://staging-ccm.phenotips.genomics4rd.ca/rest/variant-source-files/metadata",
            params=params,
            **self.base_request_args,
        )

        returned_record_count = res.json()["meta"]["returned"]

        content = res.json().get("data", [])

        found = [
            record
            for record in content
            if record["attributes"]["patientId"] == patient_id
        ]

        if found:
            return found[0]["attributes"]
        elif returned_record_count > 0:
            params["patientLimit"] += PAGE_SIZE
            params["patientOffset"] += PAGE_SIZE
            return self.get_latest_job_metadata_for_patient(
                patient_id, complete, params
            )
        else:
            raise Exception("NOT FOUND!")

    def get_patient_external_id_by_internal_id(self, patient_id: str) -> str:
        """get the patient's external ID from the internal ID"""
        res = requests.get(
            f"{self.base_url}/rest/patients/{patient_id}",
            **self.base_request_args,
        )
        return res.json()["external_id"]

    def get_variant_count(self, max=5000) -> str:
        """return count of all variants in the store up to max"""
        params = {"limit": max}
        res = requests.get(
            f"{self.base_url}/rest/variants", params=params, **self.base_request_args
        )
        if res.ok:
            return res.json()["meta"]["returned"]
        else:
            return res.get("status_code", "unknown error")


def rename_callset_fields(fieldName: str):
    if (
        fieldName.startswith("Zygosity.")
        or fieldName.startswith("Burden.")
        or fieldName.startswith("Alt_depths.")
    ):
        return sub("\..+", "", fieldName)
    else:
        return fieldName


def clean_report(filepath: str) -> str:
    """assumes a singleton, cleans the report, saves it with a new name in the same directory, and returns new path"""
    report = pd.read_csv(filepath)
    report.rename({"Position": "#Position"}, axis=1, inplace=True)
    report.rename(rename_callset_fields, axis=1, inplace=True)
    report.loc[
        (report["Zygosity"] != "Het") & (report["Zygosity"] != "Hom"), "Zygosity"
    ] = None
    missing_cols = set(ALLOWED_FIELDS).difference(set(report.columns.values))
    for col in missing_cols:
        report[col] = None
    extra_cols = set(report.columns.values).difference(set(ALLOWED_FIELDS))
    report.drop(columns=extra_cols, inplace=True)

    # reorder columns, in case it matters
    report = report[ALLOWED_FIELDS]

    saved_path = f"{path.splitext(filepath)[0]}-formatted.csv"

    with open(saved_path, "w") as f:
        f.write(report.to_csv(index=False))

    return saved_path
