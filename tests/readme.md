# Tests

Panels of tests to evaluate, optimize, and regression test sPCR.

## Finding test datasets

Ideally the datasets are small to make download easy. This directory includes a python script to query SRA for small datasets for a specific clade, for example:

     python sra.py --organism Hexacorallia

## Adding a new test dataset

Add it to `config.yaml`.

## Running tests

To run all tests, execute `batch_snake.sh`.