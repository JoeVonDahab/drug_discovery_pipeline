import os
import pickle
from typing import List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ast
import requests
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem.QED import qed as rdkit_qed
from rdkit.Chem.QED import qed
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import Descriptors
import re

# --- Load environment variables from .env file ---
# Ensure your .env file is in the same directory as this script
load_dotenv()

# --- Access your API key ---
API_KEY = os.getenv("API_KEY")

if API_KEY:
    print("NVIDIA API Key loaded successfully.")
else:
    print("NVIDIA API Key not found. Make sure it's set in your .env file or as a system environment variable.")
    raise ValueError("API_KEY is not set. Please create a .env file with your API_KEY.")

# --- API Setup ---
invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate"
headers = {
    "Authorization": f"Bearer {API_KEY}",
    "Accept": "application/json",
}
session = requests.Session()

def tanimoto_similarity(smiles: str, reference: str) -> float:

    # Get fingerprint params
    fingerprint_radius_param = 2
    fingerprint_nbits = 2048

    # Handle the reference molecule
    reference_mol = Chem.MolFromSmiles(reference)
    if reference_mol is None:
        print(f"Warning: Invalid reference SMILES: {reference}")
        return 0.0  # Cannot calculate similarity without a valid reference

    reference_fingerprint = GetMorganFingerprintAsBitVect(
        reference_mol, radius=fingerprint_radius_param, nBits=fingerprint_nbits
    )

    # Validate the other molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0  # Invalid SMILES, similarity is 0

    fingerprint = GetMorganFingerprintAsBitVect(mol, radius=fingerprint_radius_param, nBits=fingerprint_nbits)

    # Calculate and return the Tanimoto similarity
    return TanimotoSimilarity(fingerprint, reference_fingerprint)

def physical_properties(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Warning: Invalid SMILES for physical properties: {smiles}")
        return None, None, None, None, None, None  # invalid SMILES

    # Calculate physical properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hba = Descriptors.NumHAcceptors(mol)
    hbd = Descriptors.NumHDonors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    return mw, logp, hba, hbd, tpsa, rotatable_bonds

def clean_filename(smiles_string: str) -> str:
    """
    Cleans a SMILES string to be used as a valid filename.
    Replaces characters that are not allowed in filenames with underscores.
    Limits the length to prevent excessively long filenames.
    """
    # Replace common problematic characters with underscores
    cleaned = re.sub(r'[\\/:*?"<>| ]', '_', smiles_string)
    # Remove any remaining non-alphanumeric, non-underscore characters
    cleaned = re.sub(r'[^a-zA-Z0-9_.-]', '', cleaned)
    # Truncate to a reasonable length to avoid OS limits
    return cleaned[:50]  # Keep first 50 characters, adjust as needed

def process_smiles_with_molmim(smis: str, output_dir: str = "molmim_results"):
    """
    Processes a single SMILES string using the NVIDIA Molmim API
    and returns the results as a Pandas DataFrame.

    """
    print(f"\n--- Processing SMILES: {smis} ---")

    # List to collect all report data for this SMILES
    all_report_data_for_smiles = []

    mol = Chem.MolFromSmiles(smis)
    if mol is None:
        print(f"Error: Invalid starting SMILES string provided: {smis}. Skipping.")
        return pd.DataFrame()  # Return an empty DataFrame

    original_qed_score = rdkit_qed(mol)
    print(f"Original QED for {smis}: {original_qed_score}")

    # Create a list of minimum similarities (as defined in your notebook)
    num_min_sims = 3
    min_sims = np.linspace(0.1, 0.7, num_min_sims)

    # Loop through each minimum similarity value
    for min_sim in min_sims:
        print(f"*************** min_sim: {min_sim} ********************")

        # Create the request payload
        payload = {
          "smi": smis,
          "algorithm": "CMA-ES",
          "num_molecules": 10,
          "property_name": "QED",
          "minimize": False,
          "min_similarity": min_sim,
          "particles": 20,
          "iterations": 2,
          "scaled_radius": 1
        }

        try:
            # Send the request and get the response
            response = session.post(invoke_url, headers=headers, json=payload)
            response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)
            response_json = response.json()
            print(f"API response received for min_sim {min_sim}.")

            # Extract the generated SMILES
            # Use ast.literal_eval to safely parse the string representation of a list of dicts
            generated_molecules = ast.literal_eval(response_json.get('molecules', '[]'))
            gen_smiles_list = [i['sample'] for i in generated_molecules if 'sample' in i]
            print(f"Generated SMILES count: {len(gen_smiles_list)}")

            # Get the molecule objects out of valid SMILES
            valid_mol_list = []
            for smiles_gen in gen_smiles_list:
                mol_gen = Chem.MolFromSmiles(smiles_gen)
                if mol_gen:
                    valid_mol_list.append(mol_gen)
                else:
                    print(f"Warning: Generated SMILES is invalid: {smiles_gen}")

            # Convert to canonical SMILES & deduplicate
            canonical_smiles = set()
            for mol_gen in valid_mol_list:
                canonical_smi = Chem.MolToSmiles(mol_gen, canonical=True)
                canonical_smiles.add(canonical_smi)
            canonical_smiles_list = list(canonical_smiles)
            print(f"Canonical unique SMILES count: {len(canonical_smiles_list)}")

            # Calculate Tanimoto similarity and QED score for each valid SMILES
            for generated_smiles in canonical_smiles_list:
                tanimoto = tanimoto_similarity(generated_smiles, smis)
                mol_gen = Chem.MolFromSmiles(generated_smiles)
                if mol_gen:
                    qed_score = qed(mol_gen)
                else:
                    qed_score = None  # Should not happen if already validated in valid_mol_list

                # Calculate Physical properties
                mw, logp, hba, hbd, tpsa, rotatable_bonds = physical_properties(generated_smiles)

                all_report_data_for_smiles.append({
                    'original_smiles': smis,
                    'generated_smiles': generated_smiles,
                    'qed_score': qed_score,
                    'tanimoto_similarity': tanimoto,
                    'min_sim_param': min_sim,
                    'mw': mw,
                    'logp': logp,
                    'hba': hba,
                    'hbd': hbd,
                    'tpsa': tpsa,
                    'rotatable_bonds': rotatable_bonds
                })

        except requests.exceptions.HTTPError as e:
            print(f"HTTP error occurred for min_sim {min_sim}: {e}")
            print(f"Response content: {response.text}")
        except requests.exceptions.ConnectionError as e:
            print(f"Connection error occurred for min_sim {min_sim}: {e}")
        except requests.exceptions.Timeout as e:
            print(f"Timeout error occurred for min_sim {min_sim}: {e}")
        except requests.exceptions.RequestException as e:
            print(f"An unexpected error occurred during API request for min_sim {min_sim}: {e}")
        except Exception as e:
            print(f"An error occurred during processing for min_sim {min_sim}: {e}")

    if all_report_data_for_smiles:
        # Create a pandas DataFrame from the collected report data
        report_df = pd.DataFrame(all_report_data_for_smiles)
        return report_df
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no results

def process_smi_file(smi_filepath: str, output_dir: str = "molmim_results"):
    """
    Reads SMILES strings from a .smi file and processes each one
    using the Molmim API. It expects the SMILES string to be the first column.
    Accumulates all results into a single DataFrame and saves it to a CSV.

    """
    if not os.path.exists(smi_filepath):
        print(f"Error: .smi file not found at {smi_filepath}")
        return

    print(f"\n--- Starting batch processing from file: {smi_filepath} ---")
    processed_count = 0
    skipped_count = 0
    all_results = []  # List to accumulate DataFrames

    with open(smi_filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):  # Skip empty lines and comment lines
                continue

            # Split the line by whitespace and take the first element as SMILES
            parts = line.split()
            if parts:
                current_smiles = parts[0]
            else:
                print(f"Warning: Line {line_num} in {smi_filepath} is empty or malformed. Skipping.")
                skipped_count += 1
                continue

            try:
                result_df = process_smiles_with_molmim(current_smiles, output_dir)
                if not result_df.empty:
                    all_results.append(result_df)
                processed_count += 1
            except Exception as e:
                print(f"Failed to process SMILES '{current_smiles}' from line {line_num}: {e}")
                skipped_count += 1

    if all_results:
        # Concatenate all DataFrames into one
        final_report_df = pd.concat(all_results, ignore_index=True)

        # Clean the filename for the combined report
        safe_filename = clean_filename(os.path.splitext(os.path.basename(smi_filepath))[0])
        csv_filename = os.path.join(output_dir, f'molmim_combined_results_{safe_filename}.csv')

        # Save the combined DataFrame as a CSV file
        os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
        final_report_df.to_csv(csv_filename, index=False)
        print(f"\n--- Combined CSV report saved as '{csv_filename}' ---")
    else:
        print("\n--- No results generated. No combined CSV file created. ---")

    print(f"\n--- Batch processing complete for {smi_filepath} ---")
    print(f"Total SMILES processed: {processed_count}")
    print(f"Total SMILES skipped due to errors: {skipped_count}")

# --- Main execution block ---
if __name__ == "__main__":
    # Example usage:
    # 1. Create a .smi file (e.g., 'my_smiles.smi') in the same directory
    # 2. Call the function with the path to your .smi file
    #    Make sure your .env file with API_KEY is present.

    # change this with the file name
    process_smi_file("FILENAME.smi")

"""
    # To manually create a input file SMILES use this:
    dummy_smi_filename = "example_smiles_with_data.smi"
    with open(dummy_smi_filename, "w") as f:
        f.write("O=[N+]([O-])c1c(Nc2cc(F)cc(Cl)c2)ccc2nonc12 ZINC000095596639 EPAS1 7.05\n")
        f.write("CCO ZINC12345678 ExampleGene 5.0\n")
        f.write("C(=O)O ZINC87654321 AnotherGene 6.5\n")
        f.write("invalid_smiles_here ZINC000000000000 INVALID 0.0\n")  # Added an invalid SMILES for testing error handling
"""

    # You can also process a single SMILES directly if needed
    # process_smiles_with_molmim("C1=CC=C(C=C1)C(=O)O")