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

import json
import os

base_path = r"C:\Users\zhaol\Downloads\drug_discovery_pipeline\Docking Pipeline\results_output"
ligands_with_confidence = []

for i in range(1000):
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

    # Check for confidence > -1.5
    if any(conf > -1.5 for conf in confidences):
        ligands_with_confidence.append(i)

# Output list of ligands that passed the threshold
output_path = os.path.join(base_path, "ligands_with_high_confidence.txt")
with open(output_path, "w") as out:
    out.write("Ligands with confidence > -1.5:\n")
    for idx in ligands_with_confidence:
        out.write(f"ligand_{idx}\n")

print(f"Ligands with confidence > -1.5: {ligands_with_confidence}")
print("\n")
print(len(ligands_with_confidence))
print("\n")