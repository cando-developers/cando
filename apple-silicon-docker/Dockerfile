FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get install -y wget curl gpg git software-properties-common && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Add conda to path
ENV PATH="/opt/conda/bin:${PATH}"

# RUN conda install -c conda-forge jupyterlab nodejs=17.9.0 sidecar
RUN conda install -c conda-forge jupyterlab nodejs=22.12.0 sidecar
RUN bash -c "$(curl -fsSL https://www.thirdlaw.tech/pkg/cando.sh)" || exit 1
# RUN cando-user-install || exit 1
RUN cando-apptainer-install || exit 1

RUN git clone https://github.com/cando-developers/cando-demos.git || exit 1

RUN chmod -R 777 /cando-demos

WORKDIR /cando-demos

RUN find . -name "*.ipynb" -exec jupyter trust {} +

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''"]

