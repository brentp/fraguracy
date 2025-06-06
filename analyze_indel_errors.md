# Indel Error Rate Analysis Summary

## Overview

Successfully analyzed indel error rates from fraguracy output files using Python with polars for data processing and plotly for interactive visualization.

## Data Processing Steps

### 1. Indel Errors File Processing

- **Input**: `test_ovl-indel-errors.bed.gz` (compressed BED format)
- **Columns**: chrom, start, end, count, length, bq_bin, hp_dist
- **Processing**: Grouped by length, bq_bin, and hp_dist, then summed counts
- **Result**: 1,887 unique combinations with indel counts

### 2. Counts File Processing  

- **Input**: `test_ovl-counts.txt` (tab-separated)
- **Columns**: read12, FR, bq_bin, read_pos, context, hp_dist, total_count, error_count, err_rate_lo, err_rate_hi
- **Processing**: Grouped by bq_bin and hp_dist, then summed total_count
- **Result**: 155 unique combinations with total counts

### 3. Error Rate Calculation

- **Method**: Joined indel counts with total counts on bq_bin and hp_dist
- **Formula**: error_rate = indel_count / total_count
- **Final dataset**: 1,887 records with error rates

## Results

### Key Findings

- **Indel Lengths**: 66 different indel lengths observed (ranging from -47 to +28)
- **Error Rate Range**: From ~10^-9 to ~10^-5 (highly variable)
- **Homopolymer Distance Impact**: Clear relationship between hp_dist and error rates
- **Base Quality Impact**: Interactive filtering by BQ bins reveals quality-dependent error patterns (defaults to high-quality 37-59 range)

### Summary Statistics by Indel Length

- **Deletions** (negative lengths): Generally lower error rates
- **Insertions** (positive lengths): More variable error rates  
- **Single base changes** (±1): Most common, moderate error rates

## Visualization

- **Interactive Plot**: `indel_error_rates_by_hp_dist.html` (with hover data)
- **Static Plot**: `indel_error_rates_by_hp_dist.png`
- **X-axis**: Homopolymer distance (hp_dist)
- **Y-axis**: Indel error rate (log scale)
- **Colors**: Different indel length categories (aggregated)
- **Indel Length Aggregation**:
  - Individual values for -3, -2, -1, 1, 2, 3
  - Aggregated ">3" for all insertions > 3 bases
  - Aggregated "<-3" for all deletions > 3 bases
- **Features**:
  - Log scale for error rates
  - Interactive hover showing detailed information (indel count, total count, bq_bin)
  - Grid for easier reading
  - Legend with "Indel Length" title showing length categories
  - Modern plotly-based visualization
  - **Line connections (default)**: Points are connected by lines for each indel length category to show trends across hp_distance
  - **Optional scatter-only mode**: Use `--no-lines` to show scatter plot without line connections
  - **BQ Bin Filtering**: Interactive buttons to filter data by base quality bins (defaults to 37-59)
    - "All BQ Bins": Show all data
    - "BQ: 37-59 (default)": Show only high-quality base calls
    - Individual BQ bin buttons: "BQ: 05-19", "BQ: 20-36", etc.

## Files Created

1. `analyze_indel_errors.py` - Main analysis script
2. `indel_error_rates_by_hp_dist.html` - Interactive plot
3. `indel_error_rates_by_hp_dist.png` - Static plot  
4. `analysis_summary.md` - This summary document

## Usage

### Command Line Interface

The script now uses argparse for flexible command line usage:

```bash
# Basic usage (with lines - default)
python3 analyze_indel_errors.py indel_errors.bed.gz counts.txt

# Scatter plot only (no lines)
python3 analyze_indel_errors.py --no-lines indel_errors.bed.gz counts.txt

# With custom output prefix
python3 analyze_indel_errors.py --output-prefix my_analysis indel_errors.bed.gz counts.txt

# Combined options (scatter plot with custom output)
python3 analyze_indel_errors.py --no-lines --output-prefix scatter_analysis indel_errors.bed.gz counts.txt

# View help for all options
python3 analyze_indel_errors.py --help
```

### Options

- `--no-lines`: Show scatter plot only, without connecting lines (default: show lines)
- `--output-prefix`: Specify custom output file prefix (default: `indel_error_rates_by_hp_dist`)
- `--help`: Show help message with usage examples

### Arguments

1. `indel_errors_file`: Input indel errors file (BED format, can be gzipped)
2. `counts_file`: Input counts file (tab-separated format)

The script processes the input files and generates both interactive (HTML) and static (PNG) plots along with summary statistics.
