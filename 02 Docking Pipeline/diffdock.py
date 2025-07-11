# import sys
# import requests
# import time

# url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
# header_auth = "Bearer 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"

# def _upload_asset(input):
#     assets_url = "https://api.nvcf.nvidia.com/v2/nvcf/assets"

#     headers = {
#         "Authorization": header_auth,
#         "Content-Type": "application/json",
#         "accept": "application/json",
#     }

#     s3_headers = {
#         "x-amz-meta-nvcf-asset-description": "diffdock-file",
#         "content-type": "text/plain",
#     }

#     payload = {
#         "contentType": "text/plain", 
#         "description": "diffdock-file"
#     }

#     response = requests.post(
#         assets_url, headers=headers, json=payload, timeout=30
#     )

#     response.raise_for_status()

#     asset_url = response.json()["uploadUrl"]
#     asset_id = response.json()["assetId"]

#     response = requests.put(
#         asset_url,
#         data=input,
#         headers=s3_headers,
#         timeout=300,
#     )

#     response.raise_for_status()
#     return asset_id


# protein_id = _upload_asset(open("receptor_clean.pdb", "rb").read())
# ligand_id = _upload_asset(open("ligand.sdf", "rb").read())

# headers = {
#     "Content-Type": "application/json",
#     "NVCF-INPUT-ASSET-REFERENCES": ",".join([protein_id, ligand_id]),
#     "Authorization": header_auth
# }

# r = requests.post(url, headers=headers, json={
#     "ligand": ligand_id,
#     "ligand_file_type": "sdf",
#     "protein": protein_id,
#     "num_poses": 20,
#     "time_divisions": 20,
#     "steps": 18,
#     "save_trajectory": True,
#     "is_staged": True
# })
    
    
# with open("response_status.txt", "w", encoding="utf-8") as f:
#     f.write(str(r))

# # Save URL
# with open("request_url.txt", "w", encoding="utf-8") as f:
#     f.write(url)

# # Save response body text
# with open("response_text.txt", "w", encoding="utf-8") as f:
#     f.write(r.text)    
    


# print(r, url, r.text)

import os
import requests
import time
import random
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---- CONFIG ----
input_dir = "sdf_output"
output_dir = "results_output"
receptor_path = "receptor_clean.pdb"

url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
header_auth = "Bearer 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"

# ---- ASSET UPLOAD FUNCTION ----
def _upload_asset(input_data):
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

    for attempt in range(5):  # retry up to 5 times
        try:
            response = requests.post(assets_url, headers=headers, json=payload, timeout=30)
            response.raise_for_status()
            asset_url = response.json()["uploadUrl"]
            asset_id = response.json()["assetId"]

            response = requests.put(asset_url, data=input_data, headers=s3_headers, timeout=300)
            response.raise_for_status()

            return asset_id
        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                wait = 2 ** attempt + random.uniform(0, 1)
                print(f"[WARN] Rate limited. Retrying after {wait:.2f}s...")
                time.sleep(wait)
            else:
                raise
    raise RuntimeError("Failed to upload asset after multiple attempts")

# ---- UPLOAD PROTEIN ONCE ----
with open(receptor_path, "rb") as f:
    protein_id = _upload_asset(f.read())
print(f"Protein uploaded: {protein_id}")

# ---- PROCESS ONE LIGAND ----
def process_ligand(idx, sdf_file):
    try:
        ligand_path = os.path.join(input_dir, sdf_file)
        out_folder = os.path.join(output_dir, f"ligand_{idx}")
        os.makedirs(out_folder, exist_ok=True)

        with open(ligand_path, "rb") as f:
            ligand_id = _upload_asset(f.read())

        print(f"Ligand {sdf_file} uploaded: {ligand_id}")

        headers = {
            "Content-Type": "application/json",
            "NVCF-INPUT-ASSET-REFERENCES": f"{protein_id},{ligand_id}",
            "Authorization": header_auth
        }

        payload = {
            "ligand": ligand_id,
            "ligand_file_type": "sdf",
            "protein": protein_id,
            "num_poses": 20,
            "time_divisions": 20,
            "steps": 18,
            "save_trajectory": True,
            "is_staged": True
        }

        for attempt in range(5):  # Retry logic for rate-limited inference
            response = requests.post(url, headers=headers, json=payload)
            if response.status_code != 429:
                break
            wait = 2 ** attempt + random.uniform(0, 1)
            print(f"[WARN] Inference rate-limited. Retrying after {wait:.2f}s...")
            time.sleep(wait)

        with open(os.path.join(out_folder, "response_status.txt"), "w") as f:
            f.write(str(response))

        with open(os.path.join(out_folder, "request_url.txt"), "w") as f:
            f.write(url)

        with open(os.path.join(out_folder, "response_text.txt"), "w") as f:
            f.write(response.text)

        print(f"Completed ligand_{idx}: {response.status_code}")
    except Exception as e:
        print(f"[ERROR] Ligand {idx} failed: {e}")

# ---- MULTITHREADING EXECUTION ----
os.makedirs(output_dir, exist_ok=True)
sdf_files = [f for f in os.listdir(input_dir) if f.endswith(".sdf")]

# Reduce concurrency to avoid 429 errors
max_workers = min(3, len(sdf_files))  # Try 2–3 threads instead of 8

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(process_ligand, idx, sdf_file) for idx, sdf_file in enumerate(sdf_files)]
    for future in as_completed(futures):
        pass  # Ensures we wait for all tasks


