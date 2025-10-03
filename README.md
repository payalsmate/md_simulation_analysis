
# MD Simulation Analysis and Plotting

This project provides a **master script** for automating the analysis and visualization of **GROMACS molecular dynamics (MD) simulations**. It integrates multiple steps: data extraction, analysis, plotting, and movie creation, into a single pipeline.

---

## Project Structure
The project contains the following main scripts:
```
md_simulation_analysis/
├─ master.sh          # Main driver script (entry point)
├─ analyze_md.sh      # Performs structural and energetic analysis using GROMACS tools
├─ plot_analysis.py   # Generates plots (protein backbone, ligand interactions, RMSD, etc.)
├─ make_movie.sh      # Creates MD trajectory movies in chunks and concatenates them
├─ config.json        # Optional configuration file for paths and parameters
├─ .gitignore         # Ignore large simulation data, temp files, logs, etc.
└─ README.md          # You are here

```
---

## Requirements
- **Software**
  - [GROMACS](http://www.gromacs.org/) (tested version 2024.3)
  - Python 3 with required libraries:
    - `matplotlib`
    - `numpy`
    - `pandas`
    - `seaborn`
  - `jq` (for parsing JSON config)
  - FFmpeg (if making movies)

- **System**
  - Linux or WSL (tested on Ubuntu-based distributions)
  - Bash shell

---

## Usage

### 0. Make all the scripts executable
```commandline
chmod +x *.sh
```

### 1. Using a Config File
You can pass all arguments via a JSON config file:

```
./master.sh --config path/to/config.json
```

Example `config.json`:

```json
{
  "input_dir": "simulations/run1",
  "output_dir": "analysis/run1_results",
  "traj_name": "md_0_1.xtc",
  "steps": 5000,
  "protein_name": "protein.pdb",
  "ligand_name": "ligand.sdf",
  "protein_group_name": "Protein",
  "ligand_group_name": "LIG",
  "run_analysis": true,
  "run_plots": true,
  "run_movie": false,
  "num_movie_chunks": 4,
  "prefix_subplot_titles": true
}
```

### 2. Command-Line Arguments (override config)

```./master.sh \
  --input-dir simulations/run1 \
  --output-dir analysis/run1_results \
  --traj-name md_0_1.xtc \
  --protein-name protein.pdb \
  --ligand-name ligand.sdf \
  --all
```

### 3. Modes of Execution

* Run all steps:

  ```
  ./master.sh --all
  ```
* Only analysis:

  ```
  ./master.sh --analysis-only
  ```
* Only plots:

  ```
  ./master.sh --plots-only
  ```
* Only movie:

  ```
  ./master.sh --movie-only
  ```

---

## Command line Options

| Option                    | Description                                             |
| ------------------------- | ------------------------------------------------------- |
| `--input-dir`             | Directory containing simulation trajectory and topology |
| `--output-dir`            | Directory to save results                               |
| `--traj-name`             | Trajectory file (e.g., `md_0_1.xtc`)                    |
| `--steps`                 | Number of frames/steps for analysis                     |
| `--protein-name`          | Protein structure file                                  |
| `--ligand-name`           | Ligand file                                             |
| `--protein-group-name`    | Protein group (GROMACS index)                           |
| `--ligand-group-name`     | Ligand group (GROMACS index)                            |
| `--compare-dir`           | Reference directory for comparative plots               |
| `--num-movie-chunks`      | Number of chunks to split trajectory for movie making   |
| `--prefix-subplot-titles` | Add prefixes to subplot titles                          |
| `--all`                   | Run analysis, plots, and movie together                 |
| `--analysis-only`         | Run only analysis                                       |
| `--plots-only`            | Run only plotting                                       |
| `--movie-only`            | Run only movie creation                                 |

---

## Outputs

Depending on options, the pipeline generates:

* **Analysis results** (RMSD, RMSF, hydrogen bonds, distances, etc.)
* **Plots** (PNG/PDF figures for protein, ligand, and comparative results)
* **Trajectory movie** (MP4 format, optionally concatenated)

---

## Notes

* `--input-dir` and `--output-dir` are **mandatory**.
* Command-line options **override** config file values.
* Ensure that GROMACS environment variables (`gmx` command) are sourced before running.

---

## Example Workflow

1. Prepare MD simulation in GROMACS.
2. Save trajectory (`.xtc`) and topology (`.tpr`, `.gro`, `.pdb`).
3. Create a config JSON (optional).
4. Run:

   ```
   ./master.sh --config config.json
   ```
5. View results in `analysis/run1_results/`.

---
