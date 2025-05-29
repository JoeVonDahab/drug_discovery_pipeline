import json
import os

base_path = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output_actives"
highest_confidences = []

for i in range(36):
    ligand_folder = os.path.join(base_path, f"ligand_{i}")
    input_file = os.path.join(ligand_folder, "response_text.txt")
    output_folder = os.path.join(ligand_folder, "diffdock_actual_outcome")

    if not os.path.exists(input_file):
        print(f"Missing file in {ligand_folder}, skipping.")
        continue

    os.makedirs(output_folder, exist_ok=True)

    with open(input_file, "r") as f:
        data = json.load(f)

    # Write PDB files
    for j, pose in enumerate(data.get("trajectory", []), start=1):
        with open(os.path.join(output_folder, f"pose_{j}.pdb"), "w") as pdb_file:
            pdb_file.write(pose)

    # Write SDF files
    for j, sdf in enumerate(data.get("ligand_positions", []), start=1):
        with open(os.path.join(output_folder, f"ligand_pose_{j}.sdf"), "w") as sdf_file:
            sdf_file.write(sdf)

    # Write confidence scores
    confidences = data.get("position_confidence", [])
    with open(os.path.join(output_folder, "pose_confidences.txt"), "w") as out_file:
        out_file.write("Rank \t Pose Confidence\n\n")
        for j, conf in enumerate(confidences, start=1):
            out_file.write(f"{j} \t {conf}\n")

    if confidences:  # Check if the list is not empty
        highest_confidence = max(confidences)
        highest_confidences.append(highest_confidence)

# Output the highest confidences
output_path = os.path.join(base_path, "highest_confidences.txt")
with open(output_path, "w") as out:
    out.write("Highest confidence score per ligand:\n")
    for score in highest_confidences:
        out.write(f"{score}\n")

# Compute average
if highest_confidences:
    average_confidence = sum(highest_confidences) / len(highest_confidences)
    print(f"Average highest confidence across ligands: {average_confidence}")
else:
    print("No confidence scores found.")

