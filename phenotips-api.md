# Notes on the Phenotips REST API

PhenoTips provides a REST API for interacting with its data. Note that not all of the endpoints are active on all PhenoTips instances.

If you have a PhenoTips account (free and automated signup, not associated with a particular instance), you can access detailed documentation on the REST API [here](https://help.phenotips.com/hc/en-us/articles/360048543632-Variant-Store-Add-on-REST-API). The following documentation summarizes relevant information from that page and provides some additional examples specific to our use cases.

### Endpoints

The production endpoint for the PT G4RD instance is `https://phenotips.genomics4rd.ca/rest/`.
The staging endpoint is `https://staging-ccm.phenotips.genomics4rd.ca/rest/`.

Currently, the G4RD production PT instance does not have the variant or matching endpoints enabled.
The staging instance is behind the CHEO-RI firewall and can be queried via the proxy host.
Note that the staging endpoint uses basic auth, while the production endpoint requires an OAuth token (see below).

### Authentication

Depending on the instance of PhenoTips, you may need to authenticate with [basic auth](https://en.wikipedia.org/wiki/Basic_access_authentication) or [resource owner password flow](https://auth0.com/docs/flows/call-your-api-using-resource-owner-password-flow)

If you are using resource owner password flow, you'll need to make a token request to the auth server. The following example targets the G4RD Auth0 instance.

```bash
curl --request POST \
  --url 'https://phenotips.auth0.com/oauth/token' \
  --header 'content-type: application/x-www-form-urlencoded' \
  --data 'grant_type=http://auth0.com/oauth/grant-type/password-realm' \
  --data 'audience=https://phenotips.genomics4rd.ca/rest/' \
  --data 'scope=openid profile email' \
  --data 'realm=genomics4rd’ \
  --data 'client_id=<client_id>' \
  --data 'username=<username>' \
  --data 'password=<password>'
```

We are working on getting dedicated user credentials for our applications. Existing individual credentials may be used for testing in the meantime.

If successful, the token request will return a result similar to the below:

```json
{“access_token”:”…”, “id_token”:”…”, "scope":"openid profile email", "expires_in":86400, "token_type":"Bearer"}
```

Note that the **id_token** is what should be included in the `Authorization` header.

### Patient Endpoint

Detailed documentation for making calls to the Patients REST API is available [here](https://help.phenotips.com/hc/en-us/articles/360046289011-Patients-REST-API)

Note that many STAGER participants may not be in the G4RD PT instance. Recently a researcher reported:

> About half the participant IDs from datasets in STAGER are present in G4RD PhenoTips. Many of the older participants seem to be in the PhenomeCentral instance.

In STAGER terms, the convention for identifying participants in the G4RD PT data store is `<family_codename>_<participant_codename>`. Thus, to retrieve data for a patient with a `participant.participant_codename` value of `FOO` and a `family.family_codename` of `BAR` in STAGER, the query would look like:

```code
https://phenotips.genomics4rd.ca/rest/patients/fetch?eid=BAR_FOO
```

or

```code
https://phenotips.genomics4rd.ca/rest/patients/eid/BAR_FOO
```

Note that the **internal** PT patient identifier begins with a P and increments, e.g., `P0000002`

### **Upload a variant file to a patient record**

To insert a variants into the PT variant store and associate them with a patient, a csv file with the relevant data shoud be sent to the following endpoint `PUT /rest/variant-source-files/patients/<patientId>/files/<fileName>`.

The `Content-type` header shold be `multipart/form-data`. The following fields should be included:

- `metadata`: the metadata as a JSON string with keys: patientId, refGenome, fileName
- `fileStream`: the contents of the text file to upload (maximum 10 MB). The file must be in the same format as provided to the PhenoTips team when the variant store is installed.
- NB: Valid filenames must match the following regex: `[a-zA-Z0-9-_.]+`
- NB: The patient record ID and file name need to be specified in both the URL and the metadata string.

```bash
curl -X PUT
  --header 'Authorization: Bearer <id_token>'  \
  --header 'Content-Type: multipart/form-data'  \
  -F 'metadata={"patientId":"P0000037","refGenome":"GRCh37","fileName":"sample-singleton-report.csv"}' \
  -F 'fileStream=@/path/to/sample-singleton-report.csv' \
  'http://staging-ccm.phenotips.genomics4rd.ca/rest/variant-source-files/patients/P0000037/files/sample-singleton-report.csv'
```

This endpoint will return a `200` if the the request was allowed into the processing job queue. However, this does not mean that the processing job was successful. You can check the status of a job by sending a GET request to the `/rest/variant-source-files/metadata` endpoint. Note that this endpoint returns statuses for all jobs, so you may need to page through results until you find the participant you are looking for.

```bash
curl --header 'authorization: Basic <credentials>'  'https://staging-ccm.phenotips.genomics4rd.ca/rest/variant-source-files/metadata?patientOffset=25&patientLimit=25'
```

### **Delete a file**

If a processing job does not succeed, the file will remain associated with the patient, and attempts to repost it will be rejected with a `409`. In order to post a new file, you will need to manually delete the old file by sending a request to the following endpoint:

`DELETE /rest/variant-source-files/patients/<patientId>/files/<fileName>`

### **Variant endpoint**

You can query variants from the store using the following endpoint and parameters.

`GET /rest/variants`

```bash
url -v -u '<username>:<password>' \
  --get \
  -d 'offset=0' \
  -d 'limit=10' \
  -d 'sort=chrom::asc' \
  -d 'sort=pos::asc' \
  -d 'filter=patient_ids::=::P0000037' \
  -d 'filter=gene::=::ABCB5||EPX' \
  'http://staging-ccm.phenotips.genomics4rd.ca/variants'

```

### **Matching endpoint**

This is a sample request to the matcing endpoint that passes internal validation. More details will be provided once sample responses are available.

```bash
curl -X POST 'https://staging-ccm.phenotips.genomics4rd.ca/rest/variants/match' --header 'Accept: application/json' \
--header 'Authorization: <authheader> \
--header 'Content-Type: application/json' \
--data-raw '{
    "gene": {
        "geneName": "SASS6",
      	"ensemblID": "ENSG00000156876"
    },
    "variant": {
        "maxFrequency": 0.05,
        "assemblyId": "GRCh37"
    }
}'
```
