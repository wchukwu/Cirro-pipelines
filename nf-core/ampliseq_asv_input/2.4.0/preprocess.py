#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset

if __name__ == "__main__":

    ds = PreprocessDataset.from_running()

    # log
    ds.logger.info(ds.params)
