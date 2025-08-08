#!/bin/bash

#Author: Payal S. Mate

#Master script for MD Simulation Analysis and Plotting

#support function to resolve relative paths relative to config directory
resolve_relative_to_config() {
  local path="$1"
  # Leave empty paths alone
  [ -z "$path" ] && echo ""
  # if absolute pat, just return it
  if [[ "$path" = /* ]]; then
    echo "$path"
  else
    # if relative path, resolve relative to config folder
    echo "$CONFIG_DIR/$path"
  fi
}

get_config_value() {
  local key="$1"
  local value
  value=$(jq -r "if has(\"$key\") and .[\"$key\"] != null then .[\"$key\"] else empty end" "$CONFIG_FILE")
  echo "$value"
}

normalize_bool() {
  local val=$(echo "$1" | tr '[:upper:]' '[:lower:]' | xargs)  # lowercase + trim
  [[ "$val" == "true" ]] && echo "true" || echo "false"
}


# Default config values
CONFIG_FILE=""
NUM_MOVIE_CHUNKS=1
RUN_ANALYSIS=false
RUN_PLOTS=false
RUN_MOVIE=false
STEPS=1
PREFIX_SUBPLOT_TITLES=false

#extract --config argument first, if present
args=("$@")
for (( i=0; i<$#; i++ )); do
  if [[ "${args[$i]}" == "--config" ]]; then
    CONFIG_FILE="${args[$((i+1))]}"
    break
  fi
done

# Load config file if provided 
if [ -n "$CONFIG_FILE" ]; then
  if [ ! -f "$CONFIG_FILE" ]; then
    echo "Config file not found: $CONFIG_FILE"
    exit 1
  fi

  echo "Loading config values from $CONFIG_FILE"

  CONFIG_DIR=$(cd "$(dirname "$CONFIG_FILE")" && pwd)

  # Load config values (default to empty if missing)
  val=$(get_config_value "input_dir")
  [ -n "$val" ] && INPUTDIR=$(resolve_relative_to_config "$val")

  val=$(get_config_value "output_dir")
  [ -n "$val" ] && OUTDIR=$(resolve_relative_to_config "$val")

  val=$(get_config_value "compare_dir")
  [ -n "$val" ] && COMPARE_DIR=$(resolve_relative_to_config "$val")

  val=$(get_config_value "traj_name")
  [ -n "$val" ] && TRAJ_NAME="$val"

  val=$(get_config_value "steps")
  [ -n "$val" ] && STEPS="$val"

  val=$(get_config_value "protein_name")
  [ -n "$val" ] && PROTEIN_NAME="$val"

  val=$(get_config_value "ligand_name")
  [ -n "$val" ] && LIGAND_NAME="$val"

  val=$(get_config_value "protein_group_name")
  [ -n "$val" ] && PROTEIN_GROUP_NAME="$val"

  val=$(get_config_value "ligand_group_name")
  [ -n "$val" ] && LIGAND_GROUP_NAME="$val"

  val=$(get_config_value "num_movie_chunks")
  [ -n "$val" ] && NUM_MOVIE_CHUNKS="$val"

  val=$(get_config_value "run_analysis")
  [ -n "$val" ] && RUN_ANALYSIS=$(normalize_bool "$val")

  val=$(get_config_value "run_plots")
  [ -n "$val" ] && RUN_PLOTS=$(normalize_bool "$val")

  val=$(get_config_value "run_movie")
  [ -n "$val" ] && RUN_MOVIE=$(normalize_bool "$val")

  val=$(get_config_value "prefix_subplot_titles")
  [ -n "$val" ] && PREFIX_SUBPLOT_TITLES=$(normalize_bool "$val")
fi

#parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --config)
      shift; shift
      ;;
    --input-dir)
      INPUTDIR="$2"
      shift; shift
      ;;
    --output-dir)
      OUTDIR="$2"
      shift; shift
      ;;
    --compare-dir)
      COMPARE_DIR="$2"
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
    --protein-name)
      PROTEIN_NAME="$2"
      shift; shift
      ;;
    --ligand-name)
      LIGAND_NAME="$2"
      shift; shift
      ;;
    --protein-group-name)
      PROTEIN_GROUP_NAME="$2"
      shift; shift
      ;;
    --ligand-group-name)
      LIGAND_GROUP_NAME="$2"
      shift; shift
      ;;
    --all)
      RUN_ANALYSIS=true
      RUN_PLOTS=true
      RUN_MOVIE=true
      shift
      ;;
    --analysis-only)
      RUN_ANALYSIS=true
      shift
      ;;
    --plots-only)
      RUN_PLOTS=true
      shift
      ;;
    --movie-only)
      RUN_MOVIE=true
      shift
      ;;
    --num-movie-chunks)
      NUM_MOVIE_CHUNKS="$2"
      shift; shift
      ;;
    --prefix-subplot-titles)
      PREFIX_SUBPLOT_TITLES=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

STEPS="${STEPS:-1}"

#Check required args
if [ -z "$INPUTDIR" ] || [ -z "$OUTDIR" ]; then
  echo "Both --input-dir and --output-dir are required."
  exit 1
fi

mkdir -p "$OUTDIR"

echo "Input directory:  $INPUTDIR"
echo "Output directory: $OUTDIR"

#Common trajectory args
TRAJ_ARGS=""
if [ -n "$TRAJ_NAME" ]; then
  TRAJ_ARGS="$TRAJ_ARGS --traj-name $TRAJ_NAME"
fi
if [ -n "$STEPS" ]; then
  TRAJ_ARGS="$TRAJ_ARGS --steps $STEPS"
fi

#Common analysis arguments
ANALYSIS_ARGS=""
if [ -n "$PROTEIN_GROUP_NAME" ]; then
  ANALYSIS_ARGS="$ANALYSIS_ARGS --protein-group-name $PROTEIN_GROUP_NAME"
fi
if [ -n "$LIGAND_GROUP_NAME" ]; then
  ANALYSIS_ARGS="$ANALYSIS_ARGS --ligand-group-name $LIGAND_GROUP_NAME"
fi
ANALYSIS_ARGS="$ANALYSIS_ARGS --batch"

#Common plot args
PLOT_ARGS=""
if [ -n "$PROTEIN_NAME" ]; then
  PLOT_ARGS="$PLOT_ARGS --protein-name $PROTEIN_NAME --plot-backbone"
fi
if [ -n "$LIGAND_NAME" ]; then
  PLOT_ARGS="$PLOT_ARGS --ligand-name $LIGAND_NAME --plot-ligand"
fi

if [ -n "$COMPARE_DIR" ]; then
  PLOT_ARGS="$PLOT_ARGS --compare-dir $COMPARE_DIR"
fi

echo "💡 PREFIX_SUBPLOT_TITLES raw: '$PREFIX_SUBPLOT_TITLES'"
if [ "$PREFIX_SUBPLOT_TITLES" = true ]; then
  PLOT_ARGS="$PLOT_ARGS --prefix-subplot-titles"
fi

#Run selected steps based on inputs

if [ "$RUN_ANALYSIS" = true ]; then
  echo "Running analysis (analyze_md.sh) "
  ANALYSIS_CMD="bash analyze_md.sh --input-dir \"$INPUTDIR\" --output-dir \"$OUTDIR\" $TRAJ_ARGS $ANALYSIS_ARGS"
  eval "$ANALYSIS_CMD"
else
  echo "Skipping analysis step"
fi

echo "$PLOT_ARGS"
if [ "$RUN_PLOTS" = true ]; then
  echo "Running plot generation (plot_analysis.py) "
  echo "Plot args: $PLOT_ARGS"
  python3 plot_analysis.py --input-dir "$OUTDIR" --output-dir "$OUTDIR" $PLOT_ARGS
else
  echo "Skipping plot generation step"
fi

if [ "$RUN_MOVIE" = true ]; then
  echo "Running movie creation (make_movie.sh) "
  MOVIE_CMD="bash make_movie.sh --input-dir \"$INPUTDIR\" --output-dir \"$OUTDIR\" $TRAJ_ARGS --num-movie-chunks $NUM_MOVIE_CHUNKS"
  if [ "$CONCAT_TRAJ" = true ]; then
    MOVIE_CMD="$MOVIE_CMD --concat"
  fi
  eval "$MOVIE_CMD"
else
  echo "Skipping movie creation step"
fi

echo "All selected steps completed."
