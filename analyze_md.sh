#!/bin/bash

#Author: Payal S. Mate

#Script to analyze molecular dynamics simulation output data from gromacs simulation

set -e

# parse input command line Arguments
BATCH_MODE=false

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --input-dir)
      INPUTDIR="$2"
      shift; shift
      ;;
    --output-dir)
      OUTDIR="$2"
      shift; shift
      ;;
    --traj-name)
      TRAJ_NAME="$2"
      shift; shift
      ;;
    --steps)
      STEPS="$2"
      shift; shift
      ;;
    --ligand-group-name)
      LIGAND_GROUP_NAME="$2"
      shift; shift
      ;;
    --protein-group-name)
      PROTEIN_GROUP_NAME="$2"
      shift; shift
      ;;
    --batch)
      BATCH_MODE=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

#set defaults
STEPS="${STEPS:-1}"
LIGAND_GROUP_NAME="${LIGAND_GROUP_NAME:-LIG}"
PROTEIN_GROUP_NAME="${PROTEIN_GROUP_NAME:-Protein}"
BACKBONE_NAME="Backbone"

#check required inputs
if [ -z "$INPUTDIR" ] || [ -z "$OUTDIR" ] || [ -z "$TRAJ_NAME" ]; then
  echo "Usage: $0 --input-dir <dir> --output-dir <dir> --traj-name <name> [--steps <N>] [--ligand-name <name>] [--protein-name <name>] [--batch]"
  exit 1
fi

# setup
mkdir -p "$OUTDIR"
echo "Input:   $INPUTDIR"
echo "Output:  $OUTDIR"
echo "Steps:   $STEPS"
echo "Traj:    $TRAJ_NAME"
echo "Protein: $PROTEIN_GROUP_NAME"
echo "Ligand:  $LIGAND_GROUP_NAME"

# Trajectory Handling
if (( STEPS == 1 )); then
  TRAJ_JOINED="$INPUTDIR/${TRAJ_NAME}.xtc"
  TPR="$INPUTDIR/${TRAJ_NAME}.tpr"
  FINAL_STRUCT="$INPUTDIR/${TRAJ_NAME}.gro"
else
  TRAJ_JOINED="$OUTDIR/${TRAJ_NAME}_full.xtc"
  TPR="$INPUTDIR/${TRAJ_NAME}_1.tpr"
  FINAL_STRUCT="$INPUTDIR/${TRAJ_NAME}_${STEPS}.gro"

  echo "Concatenating $STEPS trajectory chunks"
  rm -f "$OUTDIR/traj_files.txt"
  for i in $(seq 1 "$STEPS"); do
    echo "$INPUTDIR/${TRAJ_NAME}_${i}.xtc" >> "$OUTDIR/traj_files.txt"
  done
  yes c | head -n "$STEPS" > "$OUTDIR/confirm_c.txt"
  gmx trjcat -f $(cat "$OUTDIR/traj_files.txt") -o "$TRAJ_JOINED" -cat -settime < "$OUTDIR/confirm_c.txt"
fi

#Center and remove PBC
echo "Centering and removing PBC"
gmx trjconv -s "$TPR" -f "$TRAJ_JOINED" -o "$OUTDIR/traj_centered.xtc" -pbc mol -center -ur compact <<< "1 0"

# Create index
INDEX="$OUTDIR/index.ndx"
gmx make_ndx -f "$FINAL_STRUCT" -o "$INDEX" << EOF
q
EOF

#Group Validation
check_group() {
  local group="$1"
  gmx select -f "$FINAL_STRUCT" -s "$TPR" -n "$INDEX" -select "group \"$group\"" &> /dev/null
}

do_backbone="y"
do_ligand="y"

if ! check_group "$PROTEIN_GROUP_NAME"; then
  echo "Protein group '$PROTEIN_GROUP_NAME' not found in index. Skipping protein analysis."
  do_backbone="n"
fi

if ! check_group "$LIGAND_GROUP_NAME"; then
  echo "Ligand group '$LIGAND_GROUP_NAME' not found in index. Skipping ligand-specific analysis."
  do_ligand="n"
fi

#Prompt if not batch
if ! $BATCH_MODE; then
  read -p "Run backbone RMSD and RMSF? [y/n] " user_backbone
  read -p "Run ligand-specific analysis (RMSD, H-bonds, distance)? [y/n] " user_ligand
  [[ $user_backbone == "n" ]] && do_backbone="n"
  [[ $user_ligand == "n" ]] && do_ligand="n"
fi


### Analysis

if [[ $do_backbone == "y" ]]; then
  echo "Compute RMSD (backbone):"
  gmx rms -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -o "$OUTDIR/rmsd.xvg" -tu ns <<< "$BACKBONE_NAME $BACKBONE_NAME"
  grep -v '^[@#]' "$OUTDIR/rmsd.xvg" > "$OUTDIR/rmsd.dat"

  echo "Compute RMSF:"
  gmx rmsf -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -o "$OUTDIR/rmsf.xvg" -res <<< "$BACKBONE_NAME"
  grep -v '^[@#]' "$OUTDIR/rmsf.xvg" > "$OUTDIR/rmsf.dat"

  echo "Compute radius of gyration:"
  gmx gyrate -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -o "$OUTDIR/gyrate.xvg" <<< "$PROTEIN_GROUP_NAME"
  grep -v '^[@#]' "$OUTDIR/gyrate.xvg" > "$OUTDIR/gyrate.dat"

  echo "Compute SASA (Protein):"
  gmx sasa -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -o "$OUTDIR/sasa.xvg" -surface 1 -output 1 <<< "$PROTEIN_GROUP_NAME $PROTEIN_GROUP_NAME"
  grep -v '^[@#]' "$OUTDIR/sasa.xvg" > "$OUTDIR/sasa.dat"
fi

if [[ $do_ligand == "y" ]]; then
  echo "Compute Ligand RMSD:"
  gmx rms -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -n "$INDEX" -o "$OUTDIR/ligand_rmsd.xvg" <<< "$LIGAND_GROUP_NAME $LIGAND_GROUP_NAME"
  grep -v '^[@#]' "$OUTDIR/ligand_rmsd.xvg" > "$OUTDIR/ligand_rmsd.dat"

  echo "Compute Ligand RMSD (fit to protein):"
  gmx rms -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -n "$INDEX" -o "$OUTDIR/ligand_rmsd_fit.xvg" <<< "$PROTEIN_GROUP_NAME $LIGAND_GROUP_NAME"
  grep -v '^[@#]' "$OUTDIR/ligand_rmsd_fit.xvg" > "$OUTDIR/ligand_rmsd_fit.dat"

  echo "Compute Ligand-Protein COM Distance"
  gmx distance -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -n "$INDEX" \
    -select "com of group $LIGAND_GROUP_NAME plus com of group $PROTEIN_GROUP_NAME" \
    -oall "$OUTDIR/ligand_protein_com_dist.xvg"
  grep -v '^[@#]' "$OUTDIR/ligand_protein_com_dist.xvg" > "$OUTDIR/ligand_protein_com_dist.dat"

  echo "Compute H-bonds (Ligand ↔ Protein):"
  printf "%s\n%s\n" "$LIGAND_GROUP_NAME" "$PROTEIN_GROUP_NAME" | \
  gmx hbond -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -n "$INDEX" -num "$OUTDIR/hbond_lig_protein.xvg"
  grep -v '^[@#]' "$OUTDIR/hbond_lig_protein.xvg" > "$OUTDIR/hbond_lig_protein.dat"

  echo "Compute Contacts (mindist):"
  gmx mindist -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -n "$INDEX" -od "$OUTDIR/contacts_lig_protein.xvg" -group <<< "$LIGAND_GROUP_NAME $PROTEIN_GROUP_NAME"
  grep -v '^[@#]' "$OUTDIR/contacts_lig_protein.xvg" > "$OUTDIR/contacts_lig_protein.dat"

  echo "Compute Ligand SASA:"
  gmx sasa -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -n "$INDEX" -o "$OUTDIR/sasa_ligand.xvg" \
    -surface "$LIGAND_GROUP_NAME" -output "$LIGAND_GROUP_NAME" <<< "$LIGAND_GROUP_NAME"
  grep -v '^[@#]' "$OUTDIR/sasa_ligand.xvg" > "$OUTDIR/sasa_ligand.dat"


  echo "Compute Hydrophobic Contacts (Ligand–Protein):"
  # Define a group of hydrophobic residues
  HYDROPHOBIC_SEL="resname ALA VAL LEU ILE MET PHE PRO TRP"

  # Create selection and output hydrophobic index file
  gmx select -s "$TPR" -f "$OUTDIR/traj_centered.xtc" -n "$INDEX" \
    -select "($HYDROPHOBIC_SEL) and group \"$PROTEIN_GROUP_NAME\" and within 0.4 of group \"$LIGAND_GROUP_NAME\"" \
    -on "$OUTDIR/hydrophobic_contacts_lig_protein.ndx" \
    -os "$OUTDIR/hydrophobic_contacts_lig_protein.xvg"

  # Optional: Track number of such atoms over time
  grep -v '^[@#]' "$OUTDIR/hydrophobic_contacts_lig_protein.xvg" > "$OUTDIR/hydrophobic_contacts_lig_protein.dat"

  echo "Compute Ligand Orientation:"

  # Define the two atoms in the ligand to use for orientation
  ORIENT_ATOM1="C1"   # Replace with actual atom name
  ORIENT_ATOM2="C5"   # Replace with actual atom name

  # Run Python orientation script
  python3 run_ligand_orientation_analysis.py \
    --traj "$OUTDIR/traj_centered.xtc" \
    --top "$FINAL_STRUCT" \
    --ligand "$LIGAND_GROUP_NAME" \
    --atom1 "$ORIENT_ATOM1" \
    --atom2 "$ORIENT_ATOM2" \
    --out "$OUTDIR/ligand_orientation.dat"
fi

echo "Analysis complete. Results saved in: $OUTDIR"
