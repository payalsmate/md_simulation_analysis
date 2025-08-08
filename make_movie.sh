#!/bin/bash
#Author: Payal S. Mate

#Script to make a movie from the trajectory of molecular dynamics simulation

set -e

#default values 
STEPS=1
CONCAT=false
NUM_MOVIE_CHUNKS=1

#parse Arguments 
while [[ $# -gt 0 ]]; do
  case $1 in
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
    --concat)
      CONCAT=true
      shift
      ;;
    --num-movie-chunks)
      NUM_MOVIE_CHUNKS="$2"
      shift; shift
      ;;
    *)
      echo "❌ Unknown option: $1"
      exit 1
      ;;
  esac
done

#Validate required inputs 
if [ -z "$INPUTDIR" ] || [ -z "$OUTDIR" ] || [ -z "$TRAJ_NAME" ]; then
  echo "Usage: ./make_movie.sh --input-dir <sim_dir> --output-dir <out_dir> --traj-name <name> [--steps N] [--concat]"
  exit 1
fi

mkdir -p "$OUTDIR"

#Define TPR and TRAJ paths 
if (( STEPS == 1 )); then
  TPR="$INPUTDIR/${TRAJ_NAME}.tpr"
  TRAJ="$INPUTDIR/${TRAJ_NAME}.xtc"
else
  TPR="$INPUTDIR/${TRAJ_NAME}_1.tpr"
  TRAJ="$OUTDIR/${TRAJ_NAME}_full.xtc"

  if $CONCAT; then
    echo "Concatenating $STEPS chunks into $TRAJ"
    rm -f "$OUTDIR/traj_files.txt"
    for i in $(seq 1 "$STEPS"); do
      echo "$INPUTDIR/${TRAJ_NAME}_${i}.xtc" >> "$OUTDIR/traj_files.txt"
    done
    yes c | head -n "$STEPS" > "$OUTDIR/confirm_c.txt"
    gmx trjcat -f $(cat "$OUTDIR/traj_files.txt") -o "$TRAJ" -cat -settime < "$OUTDIR/confirm_c.txt"
  else
    echo "Skipping concatenation — using existing: $TRAJ"
  fi
fi

#  Generate movie-ready trajectory 
echo "Checking trajectory duration with gmx check"

CHECK_OUTPUT=$(gmx check -f "$TRAJ" 2>&1)
echo "$CHECK_OUTPUT" > "$OUTDIR/debug_gmx_check.txt"

# Extract time from 'Last frame' line using regex-safe approach
TOTAL_TIME_PS=$(echo "$CHECK_OUTPUT" | awk '/Last frame/ {print int($NF)}')

if [[ -z "$TOTAL_TIME_PS" || "$TOTAL_TIME_PS" -le 0 ]]; then
  echo "Could not determine total time from gmx check output (found: '$TOTAL_TIME_PS')"
  exit 1
fi

TOTAL_TIME_NS=$(( TOTAL_TIME_PS / 1000 ))
echo "Total simulation time: ${TOTAL_TIME_NS} ns"

# Determine chunk size
CHUNK_DURATION_NS=$(( TOTAL_TIME_NS / NUM_MOVIE_CHUNKS ))
if [[ "$CHUNK_DURATION_NS" -le 0 ]]; then
  echo "Computed chunk duration is zero. Check trajectory and chunk count."
  exit 1
fi


echo "Splitting movie into $NUM_MOVIE_CHUNKS chunks (~${CHUNK_DURATION_NS} ns each)"

if [[ "$NUM_MOVIE_CHUNKS" -eq 1 ]]; then
  echo "Creating single movie.pdb (full trajectory)"
  gmx trjconv -s "$TPR" -f "$TRAJ" -o "$OUTDIR/movie.pdb" -pbc mol -center -ur compact <<< "1 0"
  echo "Movie file saved at: $OUTDIR/movie.pdb"
else
  CHUNK_DURATION_NS=$(( TOTAL_TIME_NS / NUM_MOVIE_CHUNKS ))
  echo "Splitting movie into $NUM_MOVIE_CHUNKS chunks (~${CHUNK_DURATION_NS} ns each)"

  for i in $(seq 0 $((NUM_MOVIE_CHUNKS - 1))); do
    start=$((i * CHUNK_DURATION_NS))
    end=$(( (i + 1) * CHUNK_DURATION_NS - 1 ))
    chunk_file="$OUTDIR/movie_chunk_${i}.pdb"

    echo "Creating chunk: $chunk_file from ${start} to ${end} ns"
    gmx trjconv -s "$TPR" -f "$TRAJ" -o "$chunk_file" -pbc mol -center -ur compact -b "$start" -e "$end" <<< "1 0"
  done

  echo "Movie chunks saved in: $OUTDIR/"
fi

echo "Open in VMD (Extensions->Movie Maker) to create a video or directly open in pymol."
