FROM python:3.9-slim

LABEL org.opencontainers.image.source=https://github.com/sysbiochalmers/gecko
LABEL version="0.1"
LABEL description="Docker image of SysBioChalmers/DKLcat adapted for SysBioChalmers/GECKO version 3"

COPY . .
RUN pip install --no-cache-dir -r requirements.txt torch@https://download.pytorch.org/whl/cpu/torch-1.9.1%2Bcpu-cp39-cp39-linux_x86_64.whl