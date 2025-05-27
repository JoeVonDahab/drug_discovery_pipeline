import json
import re
      
#-----------------------------------------------------------------------

with open("response_text.txt", "r") as f:
    data = json.load(f)


pdb_poses = data["trajectory"]


for i, pose in enumerate(pdb_poses, start=1):
    with open(f"diffdock_actual_outcome/pose_{i}.pdb", "w") as pdb_file:
        pdb_file.write(pose)

#---------------------------------------------------------------------------

with open("response_text.txt", "r") as f:
    data = json.load(f)


sdf_entries = data["ligand_positions"]

for i, sdf in enumerate(sdf_entries, start=1):
    with open(f"diffdock_actual_outcome/ligand_pose_{i}.sdf", "w") as sdf_file:
        sdf_file.write(sdf)

#---------------------------------------------------------------------------------

with open("response_text.txt", "r") as f:
    data = json.load(f)


confidences = data["position_confidence"]


with open("diffdock_actual_outcome/pose_confidences.txt", "w") as out_file:
    out_file.write("Rank \t Pose Confidence\n\n")
    for i, conf in enumerate(confidences, start=1):
        out_file.write(f"{i} \t {conf}\n")

