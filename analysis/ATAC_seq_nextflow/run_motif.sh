#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --job-name=motif
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --output=slurm-%j.out

PROJECT_DIR="toy_data/FHL2_ATACseq"
OUT_DIR="${PROJECT_DIR}/downstream/deseq2_consensus"

SIG_UP_CSV="${OUT_DIR}/KO_vs_WT_consensus_peaks_DESeq2.sig_up.csv"
SIG_DOWN_CSV="${OUT_DIR}/KO_vs_WT_consensus_peaks_DESeq2.sig_down.csv"

BED_UP="${OUT_DIR}/sig_up.bed"
BED_DOWN="${OUT_DIR}/sig_down.bed"

MOTIF_UP_DIR="${OUT_DIR}/motif_up"
MOTIF_DOWN_DIR="${OUT_DIR}/motif_down"

THREADS="${THREADS:-32}"
GENOME="${GENOME:-mm10}"

if ! command -v findMotifsGenome.pl >/dev/null 2>&1; then
  echo "ERROR: findMotifsGenome.pl not found. Please load/install HOMER first."
  exit 1
fi

if [[ ! -f "${SIG_UP_CSV}" || ! -f "${SIG_DOWN_CSV}" ]]; then
  echo "ERROR: Missing input CSV files:"
  echo "  ${SIG_UP_CSV}"
  echo "  ${SIG_DOWN_CSV}"
  exit 1
fi

python - "${SIG_UP_CSV}" "${BED_UP}" "${SIG_DOWN_CSV}" "${BED_DOWN}" <<'PY'
import csv
import sys

sig_up_csv, bed_up, sig_down_csv, bed_down = sys.argv[1:5]

def csv_to_bed(csv_file, bed_file):
    with open(csv_file, newline="") as f, open(bed_file, "w") as o:
        reader = csv.DictReader(f)
        for row in reader:
            o.write(f"{row['Chr']}\t{row['Start']}\t{row['End']}\t{row['PeakID']}\n")

csv_to_bed(sig_up_csv, bed_up)
csv_to_bed(sig_down_csv, bed_down)
PY

findMotifsGenome.pl "${BED_UP}" "${GENOME}" "${MOTIF_UP_DIR}" -size given -mask -p "${THREADS}"
findMotifsGenome.pl "${BED_DOWN}" "${GENOME}" "${MOTIF_DOWN_DIR}" -size given -mask -p "${THREADS}"

echo "Done."
echo "KO-up motif report:   ${MOTIF_UP_DIR}/knownResults.html"
echo "WT-up motif report:   ${MOTIF_DOWN_DIR}/knownResults.html"
