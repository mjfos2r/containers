# Start with a Python base image
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    unzip \
    ca-certificates \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Install Diamond 2.1.11
WORKDIR /tmp
RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz \
    && tar xzf diamond-linux64.tar.gz \
    && mv diamond /usr/local/bin/ \
    && chmod +x /usr/local/bin/diamond \
    && rm -f diamond-linux64.tar.gz

# Install BLAST+ 2.16.0
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz \
    && wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz.md5 \
    && md5sum -c ncbi-blast-2.16.0+-x64-linux.tar.gz.md5 || exit 1 \
    && tar xzf ncbi-blast-2.16.1+-x64-linux.tar.gz \
    && cp ncbi-blast-2.16.1+/bin/* /usr/local/bin/ \
    && rm -rf ncbi-blast-2.16.1+-x64-linux.tar.gz ncbi-blast-2.16.1+/ \

# Create a working directory
WORKDIR /app

# Install Python dependencies
RUN pip install --no-cache-dir pyyaml pandas matplotlib scikit-learn

# Clone the HGTector repository
RUN git clone https://github.com/qiyunlab/HGTector.git

# Install HGTector using setup.py
WORKDIR /app/HGTector
RUN python setup.py install

# Return to app directory
WORKDIR /app

# Set environment variables to ensure Python outputs are sent directly to terminal
ENV PYTHONUNBUFFERED=1

# Test that all of the tools installed correctly and are in path!
RUN for cmd in hgtector blastn blastp; do \
    command -v cmd &>/dev/null || exit 1; \
done

# Default command when running the container
CMD ["bash"]
