## phenotips-variant-upload.md

`report_reshape_upload.py` takes CRE reports, normalizes them (because of formatting changes over the years), demultiplexes it into participant-wise reports and POSTs these to Phenotips G4RD variant store.
The staging instance uses basic auth whereas the production instance uses auth0. Therefore the authentication method passed in will determine which instance the reports will be uploaded to. This may change in the future as the staging instance also uses auth0..

### Usage 


If `--authentication-method == 'auth0'` then `--auth-token-data` is expected to be passed in. Similarly, if `--authentication-method == 'basic'` then `--username` and `--password` is expected to be passed in.

Note: The new staging instance will also use auth0 so the script will have to be changed to accept a URL to allow the user to differentiate prod from staging. 

```
usage: report_reshape_upload.py [-h] [--username USERNAME] [--password PASSWORD] [--report-dir REPORT_DIR] [--report-path REPORT_PATH] --authentication-method AUTHENTICATION_METHOD [--auth0-token-data AUTH0_TOKEN_DATA] [--resume-from-json RESUME_FROM_JSON]
                                [--mapping-file MAPPING_FILE]

optional arguments:
  -h, --help            show this help message and exit
  --username USERNAME   basic auth username
  --password PASSWORD   basic auth password
  --report-dir REPORT_DIR
                        path to folder of reports
  --report-path REPORT_PATH
                        path to a single report
  --authentication-method AUTHENTICATION_METHOD
                        Basic or Auth0
  --auth0-token-data AUTH0_TOKEN_DATA
                        body required for auth0 authentication
  --resume-from-json RESUME_FROM_JSON
                        skips external ids found in this json
  --mapping-file MAPPING_FILE
                        use mapping dataframe instead of querying Phenotips for internal identifiers

```
Please note `--report-dir` and `--report-path` are mutually exclusive, the former expects a folder of reports whereas the latter expects a singular report. 

The following is an example of uploading reports where `--report-dir` is passed in, ie. a directory of reports is used. This would be used to populate the PT variant store retroactively. 

```
python report_reshape_upload.py \
    --authentication-method='auth0' \
    --report-dir=./2022-02-25/reports-2022-02-25 \
    --auth0-token-data='{ "client_id":"",
                  "audience": "https://phenotips.genomics4rd.ca/rest/",
                  "grant_type": ",
                  "realm": "",
                  "username" : "",
                  "password" : "" ,
                  "scope": "" }'
```

The following is an example of an upload to the production instance for a singular report using `--report-path`. The idea for this specific functionality is to integrate it with pipeline automation.

```
python report_reshape_upload.py \
    --authentication-method='auth0' \
    --report-path=./2022-02-25/reports-2022-02-25/1x/106/106.wes.2019-07-24.csv \
    --auth0-token-data='{ "client_id":"",
                  "audience": "https://phenotips.genomics4rd.ca/rest/",
                  "grant_type": ",
                  "realm": "",
                  "username" : "",
                  "password" : "" ,
                  "scope": "" }'
```

The argument `--mapping-file=<mapping.csv>` can be added to query a dataframe mapping internal identifiers with `{family_codename}_{participant_codename}`, corresponding to `external_id` and Phenotips identifiers corresponding to `report_id`.  

### Logs

Logs from python's `logging` module are written to `variant-upload-{yyyy-mm-dd}.log`.

Several attributes related to the POSTing of a demultiplexed report is saved to a json. For example, the corresponding internal identifier, extra/missing columns relative to the schema, etc. An example can be seen below. The json is normalized into a tidy dataframe for additional analyses. This is helpful for identifying which patients in the Phenotips instance had variants successfully matched to an external id as well as had a report POSTed. 

If `--report-path` is passed in, then these attributes will be saved to 

`variant-store-results-{name-of-report.csv}-on-{yyyy-mm-dd}.(csv|json)`. 

If `--report-dir` is passed in, then these attributes will be saved to 

`variant-store-results-{yyyy-mm-dd}.(csv|json)`. 

```
[
    {
        "report_name": "./test/test_prod/9005.failed.subset.wes.regular.2021-02-16.csv",
        "family": "9005",
        "participants": [
            {
                "eid": "9005_HA9000",
                "iid": "P0004384",
                "variants_found": 0,
                "missing_cols": [
                    "uce_100bp",
                    "uce_200bp"
                ],
                "extra_cols": [],
                "post_status_code": 409
            },
            .....
        ]
    },
    {
        "report_name": "./test/test_prod/9000.wes.subset.regular.2020-04-17.csv",
        "family": "9000",
        "participants": [
            {
                "eid": "9000_CH9000",
                "iid": "P0004376",
                "variants_found": 0,
                "missing_cols": [
                    "uce_100bp",
                    "spliceai_score",
                    "spliceai_impact",
                    "uce_200bp"
                ],
                "extra_cols": [
                    "splicing"
                ],
                "post_status_code": 403
            },
            {
                "eid": "9000_CH9001",
                "iid": "P0004377",
                "variants_found": 0,
                "missing_cols": [
                    "uce_100bp",
                    "spliceai_score",
                    "spliceai_impact",
                    "uce_200bp"
                ],
                "extra_cols": [
                    "splicing"
                ],
                "post_status_code": 409
            },
    ....

]

```

|report_name                                                   |family|eid        |iid     |variants_found|missing_cols                                                   |extra_cols  |post_status_code|
|--------------------------------------------------------------|------|-----------|--------|--------------|---------------------------------------------------------------|------------|----------------|
|./test/test_prod/9005.failed.subset.wes.regular.2021-02-16.csv|9005  |9005_HA9000|P0004384|0.0           |['uce_100bp', 'uce_200bp']                                     |[]          |200.0           |
|./test/test_prod/9000.wes.subset.regular.2020-04-17.csv       |9000  |9000_CH9000|P0004376|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9000.wes.subset.regular.2020-04-17.csv       |9000  |9000_CH9001|P0004377|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9000.wes.subset.regular.2020-04-17.csv       |9000  |9000_CH9002|P0004378|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9001.wes.subset.regular.2021-03-19.csv       |9001  |9001_CH9003|P0004379|0.0           |['uce_100bp', 'uce_200bp']                                     |[]          |403.0           |
|./test/test_prod/9003.failed.subset.wes.2019-09-27.csv        |9003  |9003_CH9005|P0004381|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9004.failed.subset.wes.2019-11-18.csv        |9004  |9004_CH9007|P0004383|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9004.failed.subset.wes.2019-11-18.csv        |9004  |9004_CH9006|P0004382|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9006.failed.subset.wes.regular.2020-09-10.csv|9006  |9006_SK9000|P0004385|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9002.wes.subset.regular.2020-09-01.csv       |9002  |9002_CH9004|P0004380|0.0           |['uce_100bp', 'spliceai_score', 'spliceai_impact', 'uce_200bp']|['splicing']|200.0           |
|./test/test_prod/9002.wes.subset.regular.2020-09-01.csv       |9999  |9999_CH9999|403     |              |                                                               |            |                |
|./test/test_prod/9002.wes.subset.regular.2020-09-01.csv       |9999  |9999_CH9998|403     |              |                                                               |            |                |



