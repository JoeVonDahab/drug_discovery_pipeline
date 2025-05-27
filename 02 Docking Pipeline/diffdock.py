import sys
import requests
import time

url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
header_auth = "Bearer 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"

def _upload_asset(input):
    assets_url = "https://api.nvcf.nvidia.com/v2/nvcf/assets"

    headers = {
        "Authorization": header_auth,
        "Content-Type": "application/json",
        "accept": "application/json",
    }

    s3_headers = {
        "x-amz-meta-nvcf-asset-description": "diffdock-file",
        "content-type": "text/plain",
    }

    payload = {
        "contentType": "text/plain", 
        "description": "diffdock-file"
    }

    response = requests.post(
        assets_url, headers=headers, json=payload, timeout=30
    )

    response.raise_for_status()

    asset_url = response.json()["uploadUrl"]
    asset_id = response.json()["assetId"]

    response = requests.put(
        asset_url,
        data=input,
        headers=s3_headers,
        timeout=300,
    )

    response.raise_for_status()
    return asset_id


protein_id = _upload_asset(open("receptor_clean.pdb", "rb").read())
ligand_id = _upload_asset(open("ligand.sdf", "rb").read())

headers = {
    "Content-Type": "application/json",
    "NVCF-INPUT-ASSET-REFERENCES": ",".join([protein_id, ligand_id]),
    "Authorization": header_auth
}

r = requests.post(url, headers=headers, json={
    "ligand": ligand_id,
    "ligand_file_type": "sdf",
    "protein": protein_id,
    "num_poses": 20,
    "time_divisions": 20,
    "steps": 18,
    "save_trajectory": True,
    "is_staged": True
})
    
    
with open("response_status.txt", "w", encoding="utf-8") as f:
    f.write(str(r))

# Save URL
with open("request_url.txt", "w", encoding="utf-8") as f:
    f.write(url)

# Save response body text
with open("response_text.txt", "w", encoding="utf-8") as f:
    f.write(r.text)    
    


print(r, url, r.text)