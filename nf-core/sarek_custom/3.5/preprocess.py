#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import os
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
        index=["sampleIndex", "sample", "lane", "dataset"],
        columns="read",
        values="file",
        column_prefix="fastq_"
    )
    assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

    # append metadata to file paths
    samples = ds.samplesheet.set_index("sample")
    manifest: pd.DataFrame = manifest.assign(
        **{k: manifest["sample"].apply(v.get) for k, v in samples.items()}
    )
    ordering = ['patient', 'sex', 'status', 'sample', 'lane', 'fastq_1', 'fastq_2']
    manifest = manifest.reindex(columns=ordering)

    # Overwrite the 'lane' column to provide a unique value per-row
    # This is necessary to account for datasets which merge multiple flowcells
    manifest = manifest.assign(lane=[
        str(i+1)
        for i in range(manifest.shape[0])
    ])

    # Set the default value for 'sex' to be NA
    manifest = manifest.assign(
        sex=manifest['sex'].fillna('NA')
    )

    # Transform status values "Normal" -> 0 and "Tumor" -> 1
    manifest = manifest.replace(
        to_replace=dict(
            status=dict(
                normal=0,
                Normal=0,
                tumor=1,
                Tumor=1
            )
        )
    )

    # Set the default status to 0
    manifest = manifest.assign(
        status=manifest["status"].fillna(0).apply(int)
    )

    # Set the default "patient" attribute as the sample
    manifest = manifest.assign(
        patient=manifest["patient"].fillna(manifest["sample"])
    )

    return manifest


if __name__ == "__main__":

    # Load the information for this dataset
    ds = PreprocessDataset.from_running()

    # Make the samplesheet
    manifest = make_manifest(ds)

    # Log the manifest
    for line in manifest.to_csv(index=None).split("\n"):
        ds.logger.info(line)

    # Write manifest
    manifest.to_csv("manifest.csv", index=None)

    assert ds.params.get("tools") is not None, "ERROR: You must select a variant calling tool."

    # Convert the tools paramater from a list to a comma delimited string
    ds.add_param('tools', ",".join(ds.params.get("tools")), overwrite=True)

    # If an intervals file was not selected, use --no_intervals
    if not ds.params.get("intervals"):
        ds.add_param("no_intervals", True)

    # If the user selected to save alignments in BAM format, set the flag
    if ds.params.get("alignment_format", "CRAM") == "BAM":
        ds.add_param("save_output_as_bam", True)
    ds.remove_param("alignment_format", force=True)

    # log all params
    ds.logger.info(ds.params)