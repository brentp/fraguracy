use plotly::layout::{Axis, GridPattern, Layout, LayoutGrid, Legend, RowOrder};
use plotly::{Plot, Scatter};
use polars::prelude::*;

use std::path::PathBuf;

pub fn plot(f: PathBuf) {
    let df = CsvReader::from_path(f)
        .expect("error reading csv")
        .has_header(true)
        .with_delimiter('\t' as u8)
        .finish()
        .unwrap(); //.finish().unwrap();

    let layout = Layout::new().grid(LayoutGrid::new().rows(1).columns(2));
    let mut plot = Plot::new();
    plot.set_layout(layout);

    let contexts: Vec<std::string::String> = df["context"]
        .unique()
        .expect("error getting unique contexts")
        .iter()
        .map(|s| std::string::String::from(s.get_str().unwrap()))
        .collect();

    dbg!(contexts);
    eprintln!("dataframe shape: {:?}", df.shape());
    // In [19]: df.with_columns(pl.all([pl.col("a").str.contains("sick"), ~pl.col("a").str.contains("sick of ")]).alias("match"))
    /*


    contexts.iter().map(|ctx| {
        //let mask = df.column("context")?.equal(ctx);
        df.select([col("context") == lit(*ctx)]);


        let sub = df.filter(col("context") == ctx);

        let mask = df
            .column("context")
            .unwrap()
            .contains(ctx)
            .
        df.filter(&mask).map(|subset| {
            eprintln!("{:?}", subset.shape());
        })
    });

    //let t1 = Scatter::new();

    dbg!(df.head(Some(10)));
    */
}
