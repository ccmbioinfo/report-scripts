"""
populate Phenotips G4RD staging instance with data

python3 populate-staging-instance.py \
    --username=user \
    --password=password \
    --phenotips-dump=files/phenotips_2022-01-13_20-14-DEV.json
    
python3 populate-staging-instance.py \
    --username=user \
    --password=password \
    --delete-only

NOTE: Phenotips is able to load data directly from json through their UI.

Will likely need to update for auth0 if the new staging instance uses it.
"""

from PTQuery import *
import argparse
import json
from tqdm import tqdm
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--username",
        type=str,
        required=True,
        help="basic auth username",
    )
    parser.add_argument(
        "--password", type=str, required=True, help="basic auth password"
    )
    parser.add_argument(
        "--phenotips-dump", type=str, required=False, help="path to phenotips json dump"
    )

    parser.add_argument(
        "--delete-only",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="whether to only delete patient records (ie. don't create any patients)",
    )

    return parser


def delete_existing_patients(PTQuery) -> None:
    """Deletes all existing records. Not added to the module since it could be dangerous if accidentally invoked on production."""

    all_pts = query.get_patient_info()["patientSummaries"]

    if len(all_pts) > 0:

        print(f"{len(all_pts)} existing records found")

        for record in tqdm(all_pts):
            query.delete_patient(record["id"])
    else:
        print("No existing records found")


def load_json(path_to_json: str) -> dict:

    with open(path_to_json) as f:
        data = json.load(f)

    print(f"{len(data)} records loaded from {path_to_json}")

    return data


if __name__ == "__main__":

    args = get_parser().parse_args()

    if not args.phenotips_dump and not args.delete_only:
        raise ValueError(
            "The path to the json must be specified if records are to be inserted, otherwise add the '--delete-only' argument"
        )

    BASE_REQUEST_ARGS = {
        "headers": {},
        "verify": False,
    }

    base_url = "https://dev.phenotips.genomics4rd.ca"

    username, password = (args.username, args.password)

    query = PTQuery(
        base_url=base_url,
        base_request_args=BASE_REQUEST_ARGS,
        username=username,
        password=password,
    )

    print("Deleting existing records...")
    delete_existing_patients(PTQuery)

    if not args.delete_only:
        phenotips_data = load_json(args.phenotips_dump)

        print("Inserting records from dump...")
        for record in tqdm(phenotips_data):

            query.create_patient(body=record)
