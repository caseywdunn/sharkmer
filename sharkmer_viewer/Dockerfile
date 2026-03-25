FROM continuumio/miniconda3

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  curl \
  ffmpeg \
  && rm -rf /var/lib/apt/lists/*

RUN conda create -n shark -c conda-forge -c bioconda \
  python==3.10 r-base r-minpack.lm r-argparse genomescope2

SHELL ["conda", "run", "-n", "shark", "/bin/bash", "-c"]

COPY . /opt/sharkmer_viewer
RUN pip install /opt/sharkmer_viewer \
  && cp /opt/sharkmer_viewer/genomescopemovie.sh /usr/local/bin/genomescopemovie.sh \
  && chmod +x /usr/local/bin/genomescopemovie.sh

ENV CONDA_DEFAULT_ENV=shark
ENV PATH=/opt/conda/envs/shark/bin:$PATH

WORKDIR /data
