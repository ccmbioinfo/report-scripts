## Usage

Note: this documentation is a WIP.
The staging instance uses basic auth and the production instance uses auth0. 

### Production
A bearer token must be obtained from the G4RD Auth0 instance (see `phenotips-api.md`)

```py
from PTQuery import *

# obtain bearer token given credentials

data = {
    "client_id":"",
    "audience":"",
    "grant_type":"",
    "realm":"",
    "username":"",
    "password":"",
    "scope":""
}

base_url = 'https://phenotips.genomics4rd.ca'

bearer_token = get_bearer_token(data)

query = PTQuery(base_url = base_url, 
                base_request_args = None, 
                bearer_token = bearer_token)
```

### Staging

The staging instance can only be accessed through the CHEO-RI tenancy by setting up a SOCKS proxy to the jump host. 

```py
from PTQuery import *
BASE_REQUEST_ARGS = {
    "verify": False,
}

username = '<username>'
password = '<password>'

base_url = 'https://dev.phenotips.genomics4rd.ca'

query = PTQuery(base_url = base_url, 
                        base_request_args = BASE_REQUEST_ARGS, 
                        username = username,  password = password)
````

### CRUD Operations

For more information on the endpoints, see `docs/phenotips-api.md` or the official [Phenotips REST documentation](https://help.phenotips.com/hc/en-us/articles/360046288971-REST-API-Overview). Most if not all of these endpoints operate on the patient record level.

#### Create
```
PTQuery.clean_and_post_report()
PTQuery.create_patient()
```
#### Retrieve 
```
PTQuery.get_patient_info()
PTQuery.get_job_metadata_for_patient()
PTQuery.get_patient_external_id_by_internal_id()
PTQuery.get_variant_count()
PTQuery.get_match()
PTQuery.get_variant_info()
```
#### Update 
```
```
#### Delete 
```
PTQuery.delete_patient()
PTQuery.delete_report()

```
