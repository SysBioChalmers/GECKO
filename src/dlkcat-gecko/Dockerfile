FROM python:3.9-slim

LABEL org.opencontainers.image.source=https://github.com/sysbiochalmers/gecko
LABEL version="0.1-multiarch"
LABEL description="Docker image of SysBioChalmers/DKLcat adapted for SysBioChalmers/GECKO version 3"

COPY . .
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install torch>=1.8.0+cpu --index-url https://download.pytorch.org/whl/cpu