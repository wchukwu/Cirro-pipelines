#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Filter out any index files that may have been uploaded
    ds.files = ds.files.loc[
        ds.files.apply(
            lambda r: r.get('readType', 'R') == 'R',
            axis=1
        )
    ]

    # Make a wide manifest
    manifest = ds.wide_samplesheet(
        index=["sampleIndex", "sample", "lane"],
        columns="read",
        values="file",
        column_prefix="fastq_"
    ).sort_values(
        by="sample"
    )

    # The user must have specified the experimental design with group, replicate, and control
    # The antibody sample annotation is optional
    msg = f"Samples must be annotated by group, replicate, and control"
    for cname in ['group', 'replicate', 'control']:
        assert cname in ds.samplesheet.columns.values, msg

    # All of the values in the `control` columns must match a value in `group`
    all_groups = ds.samplesheet["group"].dropna().drop_duplicates().tolist()
    all_controls = ds.samplesheet["control"].dropna().drop_duplicates().tolist()
    ds.logger.info(f"Groups: {', '.join(all_groups)}")
    ds.logger.info(f"Controls: {', '.join(all_controls)}")
    for control in all_controls:
        assert control in all_groups, f"Control '{control}' not found in group column"

    # Add those values to the manifest
    manifest = manifest.assign(**{
        cname: manifest["sample"].apply(cvals.get).fillna('')
        for cname, cvals in ds.samplesheet.set_index(
            "sample"
        ).reindex(
            columns=[
                "group",
                "replicate",
                "antibody",
                "control"
            ]
        ).items()
    })
    ds.logger.info(manifest.to_csv(index=None))
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # All replicates must be integers
    ds.logger.info("Making sure that all replicate values are integers")
    manifest = manifest.assign(
        replicate = manifest["replicate"].apply(
            lambda v: int(float(v))
        )
    )
    ds.logger.info(manifest.to_csv(index=None))

    # Only keep the expected columns, combining group + replicate -> sample
    ds.logger.info("Rearranging the columns")
    manifest = manifest.assign(
        sample=lambda d: d.apply(
            lambda r: f"{r['group']}_REP{r['replicate']}",
            axis=1
        ),
        control=lambda d: d.apply(
            lambda r: "" if pd.isnull(r['control']) or r['control'] == "" else f"{r['control']}_REP{r['replicate']}",
            axis=1
        )
    ).reindex(
        columns=["sample", "fastq_1", "fastq_2", "antibody", "control"]
    )
    ds.logger.info(manifest.to_csv(index=None))

    return manifest


if __name__ == "__main__":

    ds = PreprocessDataset.from_running()
    manifest = make_manifest(ds)

    # Save the manifest
    manifest.to_csv("design.csv", index=None)

    # Add the param for the manifest
    ds.add_param("input", "design.csv")

    # log
    ds.logger.info(ds.params)
