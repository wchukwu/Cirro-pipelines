#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Make a wide manifest
manifest = ds.wide_samplesheet(
    index=["sampleIndex", "sample", "lane"],
    columns="read",
    values="file",
    column_prefix="fastq_"
).sort_values(
    by="sample"
)
ds.logger.info(manifest.to_csv(index=None))
assert manifest.shape[0] > 0, "No files detected -- there may be an error with data ingest"

# Explicitly set the input paths in the manifest
ds.add_param(
    "input_paths",
    [
        [r['sample'], [r['fastq_1'], r['fastq_2']]]
        for _, r in manifest.iterrows()
    ]
)

# Ignore `--input` as otherwise the parameter validation will throw an error
ds.add_param("schema_ignore_params", "input")

# log
ds.logger.info(ds.params)

