#!/usr/bin/env python3
"""
Analyze indel error rates from fraguracy output files.
"""

import polars as pl
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import gzip
import argparse
from pathlib import Path

def read_indel_errors(filepath):
    """Read indel errors file and sum counts by length, bq_bin, hp_dist"""
    
    # Read the file and skip comment lines but keep header
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            lines = f.readlines()
    else:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    
    # Find header line (starts with #)
    header_line = None
    data_lines = []
    for line in lines:
        # Ensure line is a string
        if isinstance(line, bytes):
            line_str = line.decode('utf-8')
        else:
            line_str = str(line)
        
        if line_str.startswith('#'):
            header_line = line_str[1:].strip()  # Remove # and whitespace
        else:
            data_lines.append(line_str.strip())
    
    # Create a temporary file-like object with header and data
    import io
    if header_line is None:
        raise ValueError("No header line found in input file")
    csv_content = str(header_line) + '\n' + '\n'.join(data_lines)
    
    df = pl.read_csv(io.StringIO(csv_content), separator='\t', ignore_errors=True, 
                   null_values=['NA', 'N/A', ''], infer_schema_length=10000)
    
    # Filter out any rows with null values in key columns
    df = df.filter(
        pl.col('length').is_not_null() & 
        pl.col('bq_bin').is_not_null() & 
        pl.col('hp_dist').is_not_null() & 
        pl.col('count').is_not_null()
    )
    
    # Group by length, bq_bin, hp_dist and sum counts
    indel_grouped = df.group_by(['length', 'bq_bin', 'hp_dist']).agg([
        pl.col('count').sum().alias('indel_count')
    ])
    
    return indel_grouped

def read_counts_file(filepath):
    """Read counts file and sum total_count by bq_bin, hp_dist"""
    
    df = pl.read_csv(filepath, separator='\t', ignore_errors=True, 
                   null_values=['NA', 'N/A', ''], infer_schema_length=10000)
    
    # Filter out any rows with null values in key columns
    df = df.filter(
        pl.col('bq_bin').is_not_null() & 
        pl.col('hp_dist').is_not_null() & 
        pl.col('total_count').is_not_null()
    )
    
    # Group by bq_bin, hp_dist and sum total_count
    counts_grouped = df.group_by(['bq_bin', 'hp_dist']).agg([
        pl.col('total_count').sum().alias('total_count')
    ])
    
    return counts_grouped

def calculate_error_rates(indel_df, counts_df):
    """Calculate indel error rates by joining indel and total counts"""
    
    # Join the dataframes on bq_bin and hp_dist
    merged = indel_df.join(counts_df, on=['bq_bin', 'hp_dist'], how='inner')
    
    # Calculate error rate
    merged = merged.with_columns([
        (pl.col('indel_count') / pl.col('total_count')).alias('error_rate')
    ])
    
    return merged

def create_plot(df, connect_lines=False):
    """Create interactive subplot with hp_dist vs error rate and indel length vs error rate"""
    
    # Aggregate indel lengths: group lengths > 3 and < -3 into single categories
    df_plot = df.with_columns([
        pl.when(pl.col('length') > 3)
        .then(pl.lit('>3'))
        .when(pl.col('length') < -3)
        .then(pl.lit('<-3'))
        .otherwise(pl.col('length').cast(pl.Utf8))
        .alias('length_category')
    ])
    
    # Data for first subplot: aggregate by length_category, bq_bin, hp_dist
    df_hp_plot = df_plot.group_by(['length_category', 'bq_bin', 'hp_dist']).agg([
        pl.col('indel_count').sum().alias('indel_count'),
        pl.col('total_count').first().alias('total_count')  # total_count should be the same for same bq_bin/hp_dist
    ]).with_columns([
        (pl.col('indel_count') / pl.col('total_count')).alias('error_rate')
    ])
    
    # Data for second subplot: aggregate by length_category and bq_bin only (sum across all hp_dists)
    df_length_plot = df_plot.group_by(['length_category', 'bq_bin']).agg([
        pl.col('indel_count').sum().alias('indel_count'),
        pl.col('total_count').sum().alias('total_count')  # sum total_count across hp_dists
    ]).with_columns([
        (pl.col('indel_count') / pl.col('total_count')).alias('error_rate')
    ])
    
    # Get unique values for filtering
    unique_bq_bins = sorted(df_hp_plot.select('bq_bin').unique().to_numpy().flatten())
    
    # Sort length categories in logical numerical order
    def sort_length_categories(categories):
        """Sort length categories in logical order: <-3, -3, -2, -1, 1, 2, 3, >3"""
        def category_sort_key(cat):
            if cat == '<-3':
                return -1000  # Sort first
            elif cat == '>3':
                return 1000   # Sort last
            else:
                return int(cat)  # Sort numerically for individual lengths
        
        return sorted(categories, key=category_sort_key)
    
    unique_categories = sort_length_categories(df_hp_plot.select('length_category').unique().to_numpy().flatten())
    
    print(f"Available BQ bins: {unique_bq_bins}")
    print(f"Length categories: {unique_categories}")
    
    # Create subplot figure
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Error Rate by Homopolymer Distance', 'Error Rate by Indel Length'),
        vertical_spacing=0.15
    )
    
    # Color palette
    colors = px.colors.qualitative.Set1
    
    # Create traces for each combination of length_category and bq_bin
    trace_info = []
    
    # First subplot: HP Distance vs Error Rate
    for i, category in enumerate(unique_categories):
        for j, bq_bin in enumerate(unique_bq_bins):
            # Filter data for this combination
            df_subset = df_hp_plot.filter(
                (pl.col('length_category') == category) & 
                (pl.col('bq_bin') == bq_bin)
            )
            
            if df_subset.height == 0:
                continue  # Skip empty combinations
            
            # Convert to numpy arrays
            x_data = df_subset.select('hp_dist').to_numpy().flatten()
            y_data = df_subset.select('error_rate').to_numpy().flatten()
            indel_count_data = df_subset.select('indel_count').to_numpy().flatten()
            total_count_data = df_subset.select('total_count').to_numpy().flatten()
            
            # Sort by hp_dist for proper line connection
            if connect_lines and len(x_data) > 1:
                sort_idx = np.argsort(x_data)
                x_data = x_data[sort_idx]
                y_data = y_data[sort_idx]
                indel_count_data = indel_count_data[sort_idx]
                total_count_data = total_count_data[sort_idx]
                mode = 'markers+lines'
            else:
                mode = 'markers'
            
            # Determine visibility (default to 37-59 only)
            visible = True if bq_bin == '37-59' else False
            
            trace_name = f'{category} (BQ: {bq_bin})'
            
            fig.add_trace(go.Scatter(
                x=x_data,
                y=y_data,
                mode=mode,
                name=trace_name,
                visible=visible,
                legendgroup=category,  # Group traces by category for consistent legend
                marker=dict(
                    color=colors[i % len(colors)],
                    size=6,
                    opacity=0.7
                ),
                line=dict(
                    color=colors[i % len(colors)],
                    width=2
                ) if connect_lines else None,
                customdata=np.column_stack((
                    indel_count_data,
                    total_count_data,
                    [bq_bin] * len(x_data)
                )),
                hovertemplate=(
                    '<b>Indel Length:</b> ' + category + '<br>' +
                    '<b>BQ Bin:</b> ' + bq_bin + '<br>' +
                    '<b>HP Distance:</b> %{x}<br>' +
                    '<b>Error Rate:</b> %{y:.2e}<br>' +
                    '<b>Indel Count:</b> %{customdata[0]}<br>' +
                    '<b>Total Count:</b> %{customdata[1]}<br>' +
                    '<extra></extra>'
                )
            ), row=1, col=1)
            
            trace_info.append({
                'bq_bin': bq_bin,
                'category': category,
                'trace_idx': len(list(fig.data)) - 1,
                'subplot': 'hp_dist'
            })
    
    # Second subplot: Indel Length vs Error Rate  
    for i, category in enumerate(unique_categories):
        for j, bq_bin in enumerate(unique_bq_bins):
            # Filter data for this combination
            df_subset = df_length_plot.filter(
                (pl.col('length_category') == category) & 
                (pl.col('bq_bin') == bq_bin)
            )
            
            if df_subset.height == 0:
                continue  # Skip empty combinations
            
            # Convert to numpy arrays
            error_rate = df_subset.select('error_rate').to_numpy().flatten()[0]
            indel_count = df_subset.select('indel_count').to_numpy().flatten()[0]
            total_count = df_subset.select('total_count').to_numpy().flatten()[0]
            
            # Determine visibility (default to 37-59 only)
            visible = True if bq_bin == '37-59' else False
            
            # For x-axis position, convert category to numeric value for plotting
            if category == '<-3':
                x_pos = -4
            elif category == '>3':
                x_pos = 4
            else:
                x_pos = int(category)
            
            trace_name = f'{category} (BQ: {bq_bin})'
            
            fig.add_trace(go.Scatter(
                x=[x_pos],
                y=[error_rate],
                mode='markers',
                name=trace_name,  # Use same name to group in legend
                visible=visible,
                legendgroup=category,  # Group with first subplot traces
                showlegend=False,  # Don't show duplicate legend entries
                marker=dict(
                    color=colors[i % len(colors)],
                    size=8,
                    opacity=0.7
                ),
                customdata=np.array([[indel_count, total_count, bq_bin]]),
                hovertemplate=(
                    '<b>Indel Length:</b> ' + category + '<br>' +
                    '<b>BQ Bin:</b> ' + bq_bin + '<br>' +
                    '<b>Error Rate:</b> %{y:.2e}<br>' +
                    '<b>Indel Count:</b> %{customdata[0]}<br>' +
                    '<b>Total Count:</b> %{customdata[1]}<br>' +
                    '<extra></extra>'
                )
            ), row=2, col=1)
            
            trace_info.append({
                'bq_bin': bq_bin,
                'category': category,
                'trace_idx': len(list(fig.data)) - 1,
                'subplot': 'length'
            })
    
    # Create buttons for BQ bin selection
    buttons = []
    
    # Add "All" button
    all_visible = [True] * len(list(fig.data))
    buttons.append(dict(
        label="All BQ Bins",
        method="update",
        args=[{"visible": all_visible}]
    ))
    
    # Add individual BQ bin buttons
    for bq_bin in unique_bq_bins:
        visible_list = []
        for trace in trace_info:
            visible_list.append(trace['bq_bin'] == bq_bin)
        
        # Mark 37-59 as default
        label = f"BQ: {bq_bin} (default)" if bq_bin == '37-59' else f"BQ: {bq_bin}"
        
        buttons.append(dict(
            label=label,
            method="update", 
            args=[{"visible": visible_list}]
        ))
    
    # Update layout with BQ bin selector
    fig.update_layout(
        title='Indel Error Rate Analysis',
        width=1000,
        height=900,
        template='plotly_white',
        legend=dict(
            title="Indel Length",
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        ),
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=buttons,
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.01,
                xanchor="left",
                y=1.02,
                yanchor="top"
            ),
        ],
        
        #annotations=[
        #    dict(text="BQ Bin Filter:", showarrow=False,
        #         x=0.01, y=1.02, yref="paper", align="left", 
        #         font=dict(size=11, color="black"))
        #]
    )
    
    # Update subplot axes
    fig.update_xaxes(title_text="Homopolymer Distance (hp_dist)", showgrid=True, gridwidth=1, gridcolor='lightgray', row=1, col=1)
    fig.update_yaxes(title_text="Indel Error Rate", type='log', showgrid=True, gridwidth=1, gridcolor='lightgray', row=1, col=1)
    
    fig.update_xaxes(title_text="Indel Length", showgrid=True, gridwidth=1, gridcolor='lightgray', row=2, col=1)
    fig.update_yaxes(title_text="Indel Error Rate", type='log', showgrid=True, gridwidth=1, gridcolor='lightgray', row=2, col=1)
    
    # Set custom x-axis labels for second subplot
    fig.update_xaxes(
        tickvals=[-4, -3, -2, -1, 1, 2, 3, 4],
        ticktext=['<-3', '-3', '-2', '-1', '1', '2', '3', '>3'],
        row=2, col=1
    )
    
    return fig

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze indel error rates from fraguracy output files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 analyze_indel_errors.py indel_errors.bed.gz counts.txt
  python3 analyze_indel_errors.py --no-lines indel_errors.bed.gz counts.txt
        """
    )
    
    parser.add_argument(
        'indel_errors_file',
        help='Input indel errors file (BED format, can be gzipped)'
    )
    
    parser.add_argument(
        'counts_file', 
        help='Input counts file (tab-separated format)'
    )
    
    parser.add_argument(
        '--no-lines',
        action='store_false',
        dest='lines',
        help='Show scatter plot only, without connecting lines (default: show lines)'
    )
    
    parser.add_argument(
        '--output-prefix',
        default='indel_error_rates_by_hp_dist',
        help='Output file prefix (default: indel_error_rates_by_hp_dist)'
    )
    
    return parser.parse_args()

def main():
    """Main analysis function"""
    
    args = parse_arguments()
    
    print(f"Reading indel errors file: {args.indel_errors_file}")
    indel_df = read_indel_errors(args.indel_errors_file)
    print(f"Indel errors shape: {indel_df.shape}")
    print("Indel errors preview:")
    print(indel_df.head())
    
    print(f"\nReading counts file: {args.counts_file}")
    counts_df = read_counts_file(args.counts_file)
    print(f"Counts shape: {counts_df.shape}")
    print("Counts preview:")
    print(counts_df.head())
    
    print("\nCalculating error rates...")
    error_rates_df = calculate_error_rates(indel_df, counts_df)
    print(f"Error rates shape: {error_rates_df.shape}")
    print("Error rates preview:")
    print(error_rates_df.head())
    
    print(f"\nCreating plot{' with connected lines' if args.lines else ' (scatter plot only)'}...")
    fig = create_plot(error_rates_df, connect_lines=args.lines)
    
    # Save the plot as HTML for interactivity
    html_output = f"{args.output_prefix}.html"
    png_output = f"{args.output_prefix}.png"
    
    fig.write_html(html_output)
    print(f"Interactive plot saved as '{html_output}'")
    
    # Try to save PNG (requires kaleido)
    try:
        fig.write_image(png_output, width=1000, height=600)
        print(f"Static plot also saved as '{png_output}'")
    except Exception as e:
        print(f"Could not save PNG (install kaleido for PNG export): {e}")
    
    # Show summary statistics with aggregated lengths
    print("\nSummary statistics (with length aggregation):")
    summary_df = error_rates_df.with_columns([
        pl.when(pl.col('length') > 3)
        .then(pl.lit('>3'))
        .when(pl.col('length') < -3)
        .then(pl.lit('<-3'))
        .otherwise(pl.col('length').cast(pl.Utf8))
        .alias('length_category')
    ])
    
    summary = summary_df.group_by('length_category').agg([
        pl.col('error_rate').mean().alias('mean_error_rate'),
        pl.col('error_rate').std().alias('std_error_rate'),
        pl.col('indel_count').sum().alias('total_indel_count'),
        pl.col('total_count').sum().alias('total_count')
    ]).sort('length_category')
    print(summary)
    
    # Show the plot
    fig.show()

if __name__ == "__main__":
    main() 