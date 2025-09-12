FROM ubuntu:22.04

# Set noninteractive mode
ENV DEBIAN_FRONTEND=noninteractive

# Update and install basic tools
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    python3 \
    python3-pip \
    python3-dev \
    python-is-python3 \
    openjdk-11-jre-headless \
    samtools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------------------------
# Install STAR 2.7.11b
# ---------------------------
RUN mkdir -p /opt/STAR && \
    cd /opt/STAR && \
    wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    ln -s /opt/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR /usr/local/bin/STAR

# ---------------------------
# Install Salmon v1.10.1
# ---------------------------
RUN cd /opt && \
    wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
    tar -xzf salmon-1.10.0_linux_x86_64.tar.gz && \
    ln -s /opt/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

# ---------------------------
# Install FastQC
# ---------------------------
RUN cd /opt && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc

# ---------------------------
# Install Fastp
# ---------------------------
RUN cd /opt && \
    wget http://opengene.org/fastp/fastp && \
    chmod +x fastp && \
    mv fastp /usr/local/bin/fastp

# ---------------------------
# Install SeqKit
# ---------------------------
RUN cd /opt && \
    wget https://github.com/shenwei356/seqkit/releases/download/v2.6.1/seqkit_linux_amd64.tar.gz && \
    tar -xzf seqkit_linux_amd64.tar.gz && \
    chmod +x seqkit && \
    mv seqkit /usr/local/bin/seqkit

# ---------------------------
# Install Python packages
# ---------------------------
RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    jupyter \
    notebook \
    scanpy \
    anndata \
    multiqc \
    biopython

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"]


