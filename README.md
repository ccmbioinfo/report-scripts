# About

Scripts for parsing participant/dataset level information from CRE reports from WES pipelines and POSTing to Phenotip and Stager's variant store.


## 1. `transform_extract_participant_dataframes.py` 

- This script takes a single report, normalizes it as report contents and formatting change over time, and separates it into participant-wise reports, where each row in the dataframe is a variant, its annotation and the participants genotype, coverage, burden, alt depth, and zygosity with respect to the variant. 
- There will be as many output reports as there are participants. 
- The output reports will have the following naming format - the date is important for identification purposes, eg. re-analyses for the same family and participants: `{family}_{participant}_{date_of_report_generation}.csv`

Usage: 

```
python3 transform_extract_participant_dataframes.py \
    --report_path=random_reports/2x/258k/258k.wes.regular.2020-04-17.csv
``` 
The above report has 3 participants and therefore yields 3 reports. 

```
258_130201A_2020-04-17.csv
258_CH0615_2020-04-17.csv
258_CH0648_2020-04-17.csv
```

## 2. `participant_report_to_stager.py`

- Takes an output participant-wise report from `transform_extract_participant_dataframes.py`, converts it into json
   * first three k:v pairs are strings for identification purposes
   * the next k:v pair is an array of objects where each object corresponds to a variant, with its corresponding genotype object nested inside each variant
- POSTs it to `api/analyses/{analysis_id}/datasets/{dataset_id}/variants`
  * the script accepts (tentatively) a dataset_id-participant_codename mapping file generated by Stager, a user ID (?), and an analysis ID, all of which are required for the endpoint

Usage: 

For testing on Stager,

```
python3 participant_report_to_stager.py \
    --participant_report_path=2000_ACH0001_2020-04-17.csv \
    --analysis_id=1 \
    --dataset_participant_mapping=test_mapping.json \
    --user_id=1
```

Example JSON:

```

{
    "analysis": "report_name",
    "participant": "family_codename_ptp_codename",
    "family": "family_codename",
    "variant": [
        {
            "chromosome": "1",
            "position": 10521657,
            "reference_allele": "G",
            "alt_allele": "A",
            "variation": "missense_variant",
            ...., 
            "genotype": {
                "zygosity": "Het",
                "burden": 1,
                "alt_depths": "68",
                "coverage": 153,
                "genotype": "G/A"
            }
        },
        {
            "chromosome": "1",
            "position": 110754450,
            "reference_allele": "AC",
            "alt_allele": "A",
            "variation": "frameshift_variant",
            ...
            "genotype": {
                "zygosity": "Het",
                "burden": 1,
                "alt_depths": "20",
                "coverage": 41,
                "genotype": "AC/A"
            }
        }
}
```





