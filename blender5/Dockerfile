FROM nvidia/cuda:13.0.0-runtime-ubuntu24.04
COPY --from=ghcr.io/astral-sh/uv:0.9.15 /uv /uvx /bin/

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="blender5" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/blender5" \
    org.opencontainers.image.description="Blender 5 with CUDA 13.0 GPU rendering support" \
    org.opencontainers.image.version="0.0.1"

# NVIDIA container runtime settings
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility,graphics
# installed tools are available immediately
ENV UV_TOOL_BIN_DIR=/usr/local/bin
# enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    ca-certificates \
    wget \
    xz-utils \
    python3 \
    python3-pip \
    # Blender runtime deps
    libxi6 \
    libxxf86vm1 \
    libxfixes3 \
    libxrender1 \
    libgl1 \
    libxkbcommon0 \
    libsm6 \
    libgomp1 \
    # add more deps here
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Install Blender 5
ARG BLENDER_VERSION=5.0.1
RUN wget -q https://download.blender.org/release/Blender5.0/blender-${BLENDER_VERSION}-linux-x64.tar.xz \
    && tar -xf blender-${BLENDER_VERSION}-linux-x64.tar.xz -C /opt \
    && mv /opt/blender-${BLENDER_VERSION}-linux-x64 /opt/blender \
    && rm blender-${BLENDER_VERSION}-linux-x64.tar.xz \
    && ln -s /opt/blender/blender /usr/local/bin/blender

# Verify installation
RUN blender --version
RUN nvidia-smi --version || echo "nvidia-smi will work at runtime with --gpus"

# Render volumes
VOLUME ["/renders", "/projects"]
WORKDIR /renders
# copy the graph import script
COPY import_3d_graph.py /usr/local/bin/import_3d_graph.py
COPY render_gpu.py /usr/local/bin/render_gpu.py
# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]