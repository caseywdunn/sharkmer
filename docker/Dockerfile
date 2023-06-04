FROM continuumio/miniconda3

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  curl \
  ffmpeg

RUN conda create -n shark -c conda-forge -c bioconda genomescope2 \
    && echo "source activate shark" > ~/.bashrc

RUN curl --proto '=https' --tlsv1.3 https://sh.rustup.rs -sSf | sh

