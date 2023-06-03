FROM continuumio/miniconda3

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  curl

RUN conda create -n shark -c conda-forge -c bioconda python==3.10 ffmpeg genomescope2 \
    && echo "source activate my_env" > ~/.bashrc

RUN curl --proto '=https' --tlsv1.3 https://sh.rustup.rs -sSf | sh

