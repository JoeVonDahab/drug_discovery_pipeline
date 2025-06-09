import os
import json
import time
import argparse
from pathlib import Path
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from huggingface_hub import InferenceClient

# Constants
DIFFDOCK_URL = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
DIFFDOCK_TOKEN = os.environ.get("DIFFDOCK_TOKEN", "Bearer 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")
RECEPTOR_PATH = Path("02 Docking Pipeline") / "receptor_clean.pdb"

AMES_PROMPT = "Given a drug SMILES string, classify its Ames mutagenicity as (A) non-mutagenic or (B) mutagenic.\\nDrug SMILES: {smiles}\\nAnswer:"
CLINTOX_PROMPT = "Given a drug SMILES string, classify whether it is clinically toxic (A) or not (B).\\nDrug SMILES: {smiles}\\nAnswer:"
BIOAVAIL_PROMPT = "Given a drug SMILES string, predict whether its oral bioavailability is < 20% (A) or >= 20% (B).\\nDrug SMILES: {smiles}\\nAnswer:"
HF_MODEL = "google/txgemma-7b-it"
HF_TOKEN = os.environ.get("HF_TOKEN")


def smiles_to_sdf(smiles: str, out_path: Path) -> Path:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol)
    writer.close()
    return out_path


def _upload_asset(data: bytes) -> str:
    assets_url = "https://api.nvcf.nvidia.com/v2/nvcf/assets"
    headers = {
        "Authorization": DIFFDOCK_TOKEN,
        "Content-Type": "application/json",
        "accept": "application/json",
    }
    s3_headers = {
        "x-amz-meta-nvcf-asset-description": "diffdock-file",
        "content-type": "text/plain",
    }
    payload = {"contentType": "text/plain", "description": "diffdock-file"}

    resp = requests.post(assets_url, headers=headers, json=payload, timeout=30)
    resp.raise_for_status()
    asset_url = resp.json()["uploadUrl"]
    asset_id = resp.json()["assetId"]

    resp = requests.put(asset_url, data=data, headers=s3_headers, timeout=300)
    resp.raise_for_status()
    return asset_id


def run_diffdock(ligand_sdf: Path, receptor_pdb: Path = RECEPTOR_PATH) -> dict:
    with open(receptor_pdb, "rb") as f:
        protein_id = _upload_asset(f.read())
    with open(ligand_sdf, "rb") as f:
        ligand_id = _upload_asset(f.read())

    headers = {
        "Content-Type": "application/json",
        "NVCF-INPUT-ASSET-REFERENCES": f"{protein_id},{ligand_id}",
        "Authorization": DIFFDOCK_TOKEN,
    }
    payload = {
        "ligand": ligand_id,
        "ligand_file_type": "sdf",
        "protein": protein_id,
        "num_poses": 20,
        "time_divisions": 20,
        "steps": 18,
        "save_trajectory": True,
        "is_staged": True,
    }
    resp = requests.post(DIFFDOCK_URL, headers=headers, json=payload)
    resp.raise_for_status()
    return resp.json()


def query_txgemma(prompt: str) -> str:
    if HF_TOKEN is None:
        raise RuntimeError("HF_TOKEN environment variable not set")
    client = InferenceClient(model=HF_MODEL, token=HF_TOKEN)
    return client.text_generation(prompt, max_new_tokens=32).strip()


def txgemma_analysis(smiles: str) -> dict:
    results = {}
    results["AMES"] = query_txgemma(AMES_PROMPT.format(smiles=smiles))
    results["ClinTox"] = query_txgemma(CLINTOX_PROMPT.format(smiles=smiles))
    results["Bioavailability"] = query_txgemma(BIOAVAIL_PROMPT.format(smiles=smiles))
    return results


def run_pipeline(smiles: str) -> None:
    temp_dir = Path("pipeline_temp")
    temp_dir.mkdir(exist_ok=True)
    sdf_path = temp_dir / "ligand.sdf"
    smiles_to_sdf(smiles, sdf_path)
    docking_result = run_diffdock(sdf_path)
    txgemma_result = txgemma_analysis(smiles)
    out = {
        "diffdock": docking_result,
        "txgemma": txgemma_result,
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run automated docking and TxGemma analysis")
    parser.add_argument("--smiles", required=True, help="Input SMILES string")
    args = parser.parse_args()
    run_pipeline(args.smiles)
