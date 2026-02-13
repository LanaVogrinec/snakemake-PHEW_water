#!/bin/bash
set -euo pipefail

cutadapt_log=$1
original_bam=$2
filtered_bam=$3
report_file=$4

# Get total trimmed reads from cutadapt log
total_trimmed_pairs=$(grep "Pairs written (passing filters)" "$cutadapt_log" | awk '{print $5}' | tr -d ',')
if [[ -z "$total_trimmed_pairs" ]]; then
  echo "ERROR: Could not extract total trimmed pairs from $cutadapt_log"
  exit 1
fi
total_trimmed_reads=$(( total_trimmed_pairs * 2 ))

# Count mapped reads in original BAM
mapped_reads=$(samtools view -c -F 4 "$original_bam")

# Count retained mapped reads in filtered BAM
filtered_reads=$(samtools view -c -F 4 "$filtered_bam")

# Calculate percentages safely
retention_percent=$(echo "scale=2; if (${mapped_reads}==0) 0 else 100*${filtered_reads}/${mapped_reads}" | bc)
viral_percent=$(echo "scale=2; if (${total_trimmed_reads}==0) 0 else 100*${filtered_reads}/${total_trimmed_reads}" | bc)

# Write report
{
  echo "Total trimmed reads (from cutadapt log, single reads): ${total_trimmed_reads}"
  echo "Reads mapped to viral representatives (BBMap): ${mapped_reads}"
  echo "Reads retained after CoverM filtering: ${filtered_reads}"
  echo "Retention (% of mapped reads): ${retention_percent}"
  echo "Viral reads (% of total trimmed reads): ${viral_percent}"
} > "$report_file"
