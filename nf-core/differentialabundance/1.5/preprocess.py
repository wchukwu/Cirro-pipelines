import csv
from cirro.helpers.preprocess_dataset import PreprocessDataset
import json


def make_samplesheet(ds: PreprocessDataset):
    samplesheet = ds.samplesheet
    # Drop rows where sample column containing 'samtools'
    samplesheet = samplesheet[~samplesheet["sample"].str.contains("samtools")]

    variable = ds.params["variable"]
    if variable not in samplesheet:
        raise ValueError(f"Column {variable} not found in samplesheet")

    # Save to a file
    samplesheet.to_csv("samplesheet.csv", index=None)

    # Set up a workflow param pointing to that file (e.g., for nf-core/rnaseq)
    ds.logger.info(samplesheet.to_csv(index=None))


def make_contrasts(ds: PreprocessDataset):
    # Generate the contrasts.csv
    variable = ds.params["variable"]
    reference = ds.params["reference"]
    target = ds.params["target"]
    ds.remove_param("reference")
    ds.remove_param("target")
    ds.remove_param("variable")

    with open("contrasts.csv", "w") as f:
        csv.writer(f).writerow(["id", "variable", "reference", "target"])
        csv.writer(f).writerow(
            [f"{reference}_vs_{target}", variable, reference, target]
        )
    with open("contrasts.csv", "r") as f:
        ds.logger.info(f.readlines())


def set_genome(ds: PreprocessDataset):
    """
    Use the genome parameter which was selected for the input dataset.
    """

    # Get the metadata set up for this dataset, which
    # includes the params of the input dataset
    input_params = ds.metadata["inputs"][0]["params"]

    ds.add_param("genome", input_params["igenomes"]["genome"])


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    make_samplesheet(ds)
    make_contrasts(ds)
    set_genome(ds)
