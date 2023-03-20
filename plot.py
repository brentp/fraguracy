import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
import polars as pl
import numpy as np


df = pl.read_csv(sys.argv[1], sep='\t')

qual_bin = "20-39"
qual_bin = "60+"

r1 = df.filter((pl.col("read12") == "r1") & (pl.col("FR") == "f") & (
    pl.col('bq_bin') == qual_bin) & (pl.col('total_count') > 0))
r2 = df.filter((pl.col("read12") == "r2") & (pl.col("FR") == "r") & (
    pl.col('bq_bin') == qual_bin) & (pl.col('total_count') > 0))

print(r1.shape, r2.shape)

contexts = list(sorted(r1['context'].unique(), reverse=False))
print(contexts)

r1 = r1.with_columns([
    (((pl.col('error_count') + 0) / (1 + pl.col('total_count')))).alias('rate')])
r2 = r2.with_columns([
    (((pl.col('error_count') + 0) / (1 + pl.col('total_count')))).alias('rate')])

r1_rate = r1['error_count'].sum() / (r1['total_count'].sum() / 3) * 1_000_000
r2_rate = r2['error_count'].sum() / (r2['total_count'].sum() / 3) * 1_000_000

cols = plotly.colors.DEFAULT_PLOTLY_COLORS

# Create figure with secondary y-axis
fig = make_subplots(
    rows=2, subplot_titles=[f"read1(F) errors  per million read-bases: {r1_rate:.3f}", f"read2(R) errors per million read-bases: {r2_rate: 3f}"],
    vertical_spacing=0.1,
)

for i, ctx in enumerate(contexts):
    sub1 = r1.filter(pl.col('context') == ctx)
    sub2 = r2.filter(pl.col('context') == ctx)

    rate1 = sub1['error_count'].sum() / sub1['total_count'].sum() * 1_000_000
    rate2 = sub2['error_count'].sum() / sub2['total_count'].sum() * 1_000_000

    t1 = go.Scatter(name=f'{ctx}', x=np.array(sub1["read_pos"]), y=(1_000_000 * np.array(
        sub1['rate'])),
        hovertemplate="rate/Mb: %{y:.2g}  <i>errors</i>:%{text}",
        text=[f'<b>{c}</b> of {n:,}' for c,
              n in zip(sub1["error_count"], sub1['total_count'])],
        line=dict(color=cols[i]))
    t2 = go.Scatter(name='f{ctx}', x=np.array(sub2["read_pos"]), y=(1_000_000 * np.array(
        sub2['rate'])),
        hovertemplate="rate/Mb: %{y:.2g}  <i>errors</i>:%{text}",
        text=[f'<b>{c}</b> of {n:,}' for c,
              n in zip(sub2["error_count"], sub2['total_count'])],
        line=dict(color=cols[i]), showlegend=False)
    fig.add_trace(t1, row=1, col=1)
    fig.add_trace(t2, row=2, col=1)

# fig.update_layout(barmode='stack')
fig.update_layout(hovermode='x unified')
fig.update_xaxes(title_text="relative read position")
# fig.update_layout(legend_traceorder="reversed")

# fig.update_layout(title_text="error-rate along a read")

fig.update_yaxes(title_text="errors per million read bases")
# fig.update_layout(yaxis_tickformat = '%g')


fig.write_html("read.html")
