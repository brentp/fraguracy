import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
import polars as pl
import numpy as np

# Create figure with secondary y-axis
fig = make_subplots(
    rows=2, subplot_titles=["read1(F)", "read2(R)"],
    vertical_spacing=0.1,
)
# rows=2, specs=[[{"secondary_y": True}], [{"secondary_y": True}]], shared_yaxes=True, shared_xaxes=True)

df = pl.read_csv(sys.argv[1], sep='\t')

qual_bin = "60+"
qual_bin = "20-39"

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

cols = plotly.colors.DEFAULT_PLOTLY_COLORS

for i, ctx in enumerate(contexts):
    sub1 = r1.filter(pl.col('context') == ctx)
    sub2 = r2.filter(pl.col('context') == ctx)

    t1 = go.Scatter(name=ctx, x=np.array(sub1["read_pos"]), y=(1_000_000 * np.array(
        sub1['rate'])), line=dict(color=cols[i]))
    t2 = go.Scatter(name=ctx, x=np.array(sub2["read_pos"]), y=(1_000_000 * np.array(
        sub2['rate'])), line=dict(color=cols[i]), showlegend=False)
    fig.add_trace(t1, row=1, col=1)
    fig.add_trace(t2, row=2, col=1)

# fig.update_layout(barmode='stack')
fig.update_xaxes(title_text="relative read position")
# fig.update_layout(legend_traceorder="reversed")

# fig.update_layout(title_text="error-rate along a read")

fig.update_yaxes(title_text="errors per million read bases")
# fig.update_layout(yaxis_tickformat = '%g')


fig.write_html("read.html")

"""
read12	FR	bq_bin	mq_bin	read_pos	context	total_count	error_count
r1	f	0-5	0-5	0	AC	0	0
r1	f	0-5	0-5	0	AG	0	0
r1	f	0-5	0-5	0	AT	0	0
r1	f	0-5	0-5	0	CA	0	0
r1	f	0-5	0-5	0	CG	0	0
r1	f	0-5	0-5	0	CT	0	0
r1	f	0-5	05-19	0	AC	0	0
r1	f	0-5	05-19	0	AG	0	0
r1	f	0-5	05-19	0	AT	0	0
"""
