# Sharkmer

## Development

### Rust

#### Conda

To build and install the conda package locally:

    conda build .
    conda install miniconda3/envs/devel/conda-bld/osx-arm64/sharkmer_rust-0.1.0-0.tar.bz2 # exact path wil come from above command


#### Debugging

To debug in vs code:

- Go to File > add folder to workspace
- Select `sharkmer/sharkmer` folder

This opens just the rust code, and debugging works fine. 