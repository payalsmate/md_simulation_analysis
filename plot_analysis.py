#Author: Payal S. Mate

#Python script to plot various metrics generated after analyzing the trajectory data of a molecular dynamics simulation

import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

# Font size settings
FONT_TITLE = 16
FONT_LABEL = 14
FONT_LEGEND = 12
FONT_TICKS = 12

def load_xy(path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


class Plotter:
    def __init__(self, prefix_subplot_titles=False):
        self.prefix_subplot_titles = prefix_subplot_titles

    def get_title(self, title, idx, override=None):
        enabled = override if override is not None else self.prefix_subplot_titles
        labels = ['[A]', '[B]', '[C]', '[D]', '[E]', '[F]', 'G']
        return f"{labels[idx]} {title}" if enabled and idx < len(labels) else title


    def plot_backbone(self, input_dir, output_dir, compare_dir=None, protein_name=None, ligand_name=None):
        plt.style.use('seaborn-v0_8-darkgrid')
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))        

        def plot_metric(ax, filename, title, xlabel, ylabel, color, compare=False):
            x, y = load_xy(os.path.join(input_dir, filename))
            main_label = f"{protein_name or 'Protein'} + {ligand_name or 'Ligand'}"
            cmp_label = f"{protein_name or 'Protein'} only"
            ax.plot(x, y, color=color, label=main_label)
            if compare and compare_dir:
                x2, y2 = load_xy(os.path.join(compare_dir, filename))
                ax.plot(x2, y2, color='gray', linestyle='--', label=cmp_label)
                ax.legend(fontsize=FONT_LEGEND)
            ax.set_title(title, fontsize=FONT_TITLE)
            ax.set_xlabel(xlabel, fontsize=FONT_LABEL)
            ax.set_ylabel(ylabel, fontsize=FONT_LABEL)
            ax.tick_params(axis='both', labelsize=FONT_TICKS)

        # RMSD
        plot_metric(axs[0, 0], 'rmsd.dat', self.get_title('Backbone RMSD',0), 'Time (ns)', 'RMSD (nm)', 'royalblue', compare=True)

        # RMSF
        x, y = load_xy(os.path.join(input_dir, 'rmsf.dat'))
        n = len(x) // 4
        colors = ['darkorange', 'crimson', 'teal', 'purple']
        for i in range(4):
            axs[0, 1].plot(x[i*n:(i+1)*n], y[i*n:(i+1)*n], label=f'Chain {chr(65+i)}', color=colors[i])
        axs[0, 1].set_title(self.get_title('RMSF per Chain',1), fontsize=FONT_TITLE)
        axs[0, 1].set_xlabel('Residue', fontsize=FONT_LABEL)
        axs[0, 1].set_ylabel('RMSF (nm)', fontsize=FONT_LABEL)
        axs[0, 1].legend(fontsize=FONT_LEGEND)
        axs[0, 1].tick_params(axis='both', labelsize=FONT_TICKS)

        # Radius of Gyration
        plot_metric(axs[1, 0], 'gyrate.dat', self.get_title('Radius of Gyration',2), 'Time (ps)', 'Rg (nm)', 'green', compare=True)

        # SASA
        plot_metric(axs[1, 1], 'sasa.dat', self.get_title('SASA (Protein)',3), 'Time (ps)', 'SASA (nm²)', 'purple', compare=True)

        plt.tight_layout()        
        # Adjust filename if comparing
        fname = 'backbone_metrics_comparison.png' if compare_dir else 'backbone_metrics.png'
        outpath = os.path.join(output_dir, fname)    
        plt.savefig(outpath, dpi=300)
        print(f"✅ Backbone plots saved to {outpath}")

    def plot_rmsf_per_chain(self, input_dir, output_dir, compare_dir=None, protein_name=None, ligand_name=None):
        plt.style.use('seaborn-v0_8-darkgrid')
        fig, axs = plt.subplots(2, 2, figsize=(14, 10))

        x, y = load_xy(os.path.join(input_dir, 'rmsf.dat'))
        if compare_dir:
            x_cmp, y_cmp = load_xy(os.path.join(compare_dir, 'rmsf.dat'))

        n = len(x) // 4  # assumes 4 chains equally split
        colors = ['darkorange', 'crimson', 'teal', 'purple']

        main_label = f"{protein_name or 'Protein'} + {ligand_name or 'Ligand'}"
        cmp_label = f"{protein_name or 'Protein'} only"

        for i in range(4):
            ax = axs[i // 2, i % 2]
            xi = x[i*n:(i+1)*n]
            yi = y[i*n:(i+1)*n]
            ax.plot(xi, yi, label= main_label, color=colors[i])

            if compare_dir:
                yi_cmp = y_cmp[i*n:(i+1)*n]
                ax.plot(xi, yi_cmp, label=cmp_label, color='gray', linestyle='--')

            ax.set_title(self.get_title(f'RMSF - Chain {chr(65+i)}', i), fontsize=FONT_TITLE)
            ax.set_xlabel('Residue', fontsize=FONT_LABEL)
            ax.set_ylabel('RMSF (nm)',fontsize=FONT_LABEL)
            ax.legend(fontsize=FONT_LEGEND)
            ax.tick_params(axis='both', labelsize=FONT_TICKS)

        plt.tight_layout()
        fname = 'rmsf_per_chain_comparison.png' if compare_dir else 'rmsf_per_chain.png'
        outpath = os.path.join(output_dir, fname)
        plt.savefig(outpath, dpi=300)
        print(f"RMSF per-chain plot saved to {outpath}")


    def plot_ligand(self, input_dir, output_dir):
        plt.style.use('seaborn-v0_8-darkgrid')

        def plot_file(ax, filename, title, xlabel, ylabel, color):
            x, y = load_xy(os.path.join(input_dir, filename))
            ax.plot(x, y, color=color)
            ax.set_title(title, fontsize=FONT_TITLE)
            ax.set_xlabel(xlabel, fontsize=FONT_LABEL)
            ax.set_ylabel(ylabel, fontsize=FONT_LABEL)
            ax.tick_params(axis='both', labelsize=FONT_TICKS)

        # --- Figure 1: RMSD, COM Distance, H-bonds ---
        fig1, axs1 = plt.subplots(2, 2, figsize=(12, 8))
        plot_file(axs1[0, 0], 'ligand_rmsd.dat', self.get_title('Ligand RMSD (self-fit)',0), 'Time (ns)', 'RMSD (nm)', 'blue')
        plot_file(axs1[0, 1], 'ligand_rmsd_fit.dat', self.get_title('Ligand RMSD (fit to protein)',1), 'Time (ns)', 'RMSD (nm)', 'navy')
        plot_file(axs1[1, 0], 'ligand_protein_com_dist.dat', self.get_title('Ligand–Protein COM Distance',2), 'Time (ps)', 'Distance (nm)', 'darkred')
        plot_file(axs1[1,1], 'contacts_lig_protein.dat', self.get_title('Ligand–Protein Contacts',3), 'Time (ps)', 'Contact Distance (nm)', 'orange')
        fig1.tight_layout()
        outpath1 = os.path.join(output_dir, 'ligand_distances.png')
        fig1.savefig(outpath1, dpi=300)
        print(f"Ligand distance-related plots saved to {outpath1}")

        # --- Figure 2: Contacts & SASA ---
        fig2, axs2 = plt.subplots(2, 2, figsize=(12, 8))
        plot_file(axs2[0,0], 'hbond_lig_protein.dat', self.get_title('H-bonds (Ligand–Protein)',0), 'Time (ps)', 'Number of H-bonds', 'seagreen')
        plot_file(axs2[0,1], 'hydrophobic_contacts_lig_protein.dat', self.get_title('Hydrophobic Contacts (Ligand–Protein)',1), 'Time (ps)', 'No. of Contacts', 'darkgreen')
        plot_file(axs2[1,0], 'sasa_ligand.dat', self.get_title('Ligand SASA',2), 'Time (ps)', 'SASA (nm²)', 'purple')
        plot_file(axs2[1,1], 'ligand_orientation.dat', self.get_title('Ligand Orientation Angle',3), 'Time (ps)', 'Angle (°)', 'teal')  # <- NEW PLOT
        fig2.tight_layout()
        outpath2 = os.path.join(output_dir, 'ligand_contacts_sasa.png')
        fig2.savefig(outpath2, dpi=300)
        print(f"Ligand contacts/SASA plots saved to {outpath2}")


#Main Function starts

# --- CLI
parser = argparse.ArgumentParser(description='Plot GROMACS analysis results.')
parser.add_argument('--input-dir', required=True, help='Directory with input .dat files')
parser.add_argument('--output-dir', required=True, help='Directory to save plots')
parser.add_argument('--plot-backbone', action='store_true', help='Plot backbone-related plots')
parser.add_argument('--plot-ligand', action='store_true', help='Plot ligand-related plots')
parser.add_argument('--protein-name', help='Name of the protein (for labeling)')
parser.add_argument('--ligand-name', help='Name of the ligand (for labeling)')
parser.add_argument('--compare-dir', help='Optional directory for comparison with protein-only sim')
parser.add_argument('--prefix-subplot-titles', action='store_true', help='Add [A], [B], etc. to subplot titles')
args = parser.parse_args()

print(' prefix_subplot_titles: ', args.prefix_subplot_titles)
plotter = Plotter(prefix_subplot_titles=args.prefix_subplot_titles)

# --- Run Selected Plots ---
print(f"Input directory: {args.input_dir}")
print(f"Output directory: {args.output_dir}")
if args.compare_dir:
    print(f"Comparing backbone data with: {args.compare_dir}")

if not args.plot_backbone and not args.plot_ligand:
    user_choice = input("Plot backbone (b), ligand (l), or both (bl)? ").strip().lower()
    if 'b' in user_choice:
        print("Plotting backbone data ")
        plotter.plot_backbone(args.input_dir, args.output_dir, args.compare_dir, args.protein_name, args.ligand_name)
        print("Plotting RMSF per chain...")
        plotter.plot_rmsf_per_chain(args.input_dir, args.output_dir, args.compare_dir, args.protein_name, args.ligand_name)
    if 'l' in user_choice:
        print("Plotting ligand data ")
        plotter.plot_ligand(args.input_dir, args.output_dir)
else:
    if args.plot_backbone:
        print("Plotting backbone data ")
        plotter.plot_backbone(args.input_dir, args.output_dir, args.compare_dir, args.protein_name, args.ligand_name)
        print("Plotting RMSF per chain ")
        plotter.plot_rmsf_per_chain(args.input_dir, args.output_dir, args.compare_dir, args.protein_name, args.ligand_name)
    if args.plot_ligand:
        print("Plotting ligand data ")
        plotter.plot_ligand(args.input_dir, args.output_dir)

print("Plotting complete.")
