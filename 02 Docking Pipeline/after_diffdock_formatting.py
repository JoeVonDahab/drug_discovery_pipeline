# import json
# import re
      
# #-----------------------------------------------------------------------

# with open("response_text.txt", "r") as f:
#     data = json.load(f)


# pdb_poses = data["trajectory"]


# for i, pose in enumerate(pdb_poses, start=1):
#     with open(f"diffdock_actual_outcome/pose_{i}.pdb", "w") as pdb_file:
#         pdb_file.write(pose)

# #---------------------------------------------------------------------------

# with open("response_text.txt", "r") as f:
#     data = json.load(f)


# sdf_entries = data["ligand_positions"]

# for i, sdf in enumerate(sdf_entries, start=1):
#     with open(f"diffdock_actual_outcome/ligand_pose_{i}.sdf", "w") as sdf_file:
#         sdf_file.write(sdf)

# #---------------------------------------------------------------------------------

# with open("response_text.txt", "r") as f:
#     data = json.load(f)


# confidences = data["position_confidence"]


# with open("diffdock_actual_outcome/pose_confidences.txt", "w") as out_file:
#     out_file.write("Rank \t Pose Confidence\n\n")
#     for i, conf in enumerate(confidences, start=1):
#         out_file.write(f"{i} \t {conf}\n")

#Average highest confidence across ligands: 0.29975198871559566

import json
import os

base_path = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output"
ligand_confidences = []

for i in range(10000):
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

    # Save highest valid confidence
    valid_confidences = [c for c in confidences if c is not None]
    if valid_confidences:
        highest = max(valid_confidences)
        ligand_confidences.append((f"ligand_{i}", highest))
    else:
        print(f"No valid confidence values in {ligand_folder}, skipping.")

# Sort ligands by highest confidence descending
ligand_confidences.sort(key=lambda x: x[1], reverse=True)

# Select as many ligands as possible with average ≥ 0.299
selected_ligands = []
running_total = 0.0

for idx, (ligand_id, conf_score) in enumerate(ligand_confidences, start=1):
    running_total += conf_score
    avg = running_total / idx
    if avg >= 0.299:
        selected_ligands.append((ligand_id, conf_score))
    else:
        break

# Output selected ligands
output_path = os.path.join(base_path, "selected_ligands_by_confidence.txt")
with open(output_path, "w") as out:
    out.write("Selected ligands with highest confidence (avg >= 0.299):\n\n")
    out.write("Ligand ID\tConfidence Score\n")
    for ligand_id, score in selected_ligands:
        out.write(f"{ligand_id}\t{score:.4f}\n")

# Print summary
if selected_ligands:
    avg_final = sum(score for _, score in selected_ligands) / len(selected_ligands)
    print(f"Selected {len(selected_ligands)} ligands with avg confidence ≥ 0.299.")
    print(f"Final average confidence: {avg_final:.4f}")
else:
    print("No ligands selected. All had average confidence below 0.299.")
