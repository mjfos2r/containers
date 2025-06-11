FROM python:3.13-slim
COPY --from=ghcr.io/astral-sh/uv:0.7.12 /uv /uvx /bin/

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*


# Set working directory
WORKDIR /app

RUN uv venv
ENV PATH="/app/.venv/bin:$PATH"
RUN uv pip install NanoPlot==1.44.1
# Verify installation
RUN command -v NanoPlot &>/dev/null || exit 1

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]