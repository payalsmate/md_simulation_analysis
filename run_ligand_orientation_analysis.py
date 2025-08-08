#Author: Payal S. Mate

#Python script to compute ligand orientation in the trajectory data of a molecular dynamics simulation

import mdtraj as md
import numpy as np
import argparse
import os

def compute_orientation(traj_file, top_file, ligand_resname, atom1_name, atom2_name, output_file):
    traj = md.load(traj_file, top=top_file)
    ligand = traj.topology.select(f"resname {ligand_resname}")
    
    atom1_index = traj.topology.select(f"resname {ligand_resname} and name {atom1_name}")
    atom2_index = traj.topology.select(f"resname {ligand_resname} and name {atom2_name}")
    
    if len(atom1_index) != 1 or len(atom2_index) != 1:
        raise ValueError("Atom selection must return exactly one atom each.")
    
    atom1 = atom1_index[0]
    atom2 = atom2_index[0]

    vec = traj.xyz[:, atom2] - traj.xyz[:, atom1]  # shape: (n_frames, 3)
    z_axis = np.array([0, 0, 1])

    # Normalize vectors
    vec_norm = vec / np.linalg.norm(vec, axis=1)[:, None]
    dot_product = np.dot(vec_norm, z_axis)
    angles_rad = np.arccos(np.clip(dot_product, -1.0, 1.0))
    angles_deg = np.degrees(angles_rad)

    times = traj.time  # in ps
    data = np.column_stack((times, angles_deg))
    np.savetxt(output_file, data, fmt="%.2f", header="Time(ps) Angle(deg)")
    print(f"Orientation data saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ligand orientation analysis")
    parser.add_argument('--traj', required=True, help='Trajectory file (e.g., .xtc)')
    parser.add_argument('--top', required=True, help='Topology file (e.g., .pdb or .gro)')
    parser.add_argument('--ligand', required=True, help='Ligand residue name')
    parser.add_argument('--atom1', required=True, help='Ligand atom name (start of vector)')
    parser.add_argument('--atom2', required=True, help='Ligand atom name (end of vector)')
    parser.add_argument('--out', required=True, help='Output file (e.g., ligand_orientation.dat)')
    args = parser.parse_args()

    compute_orientation(args.traj, args.top, args.ligand, args.atom1, args.atom2, args.out)
