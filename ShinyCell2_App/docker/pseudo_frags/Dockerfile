FROM python:3.10-slim-bullseye
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential pkg-config libhdf5-dev libatlas-base-dev liblapack-dev gfortran \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY h5mu_to_pseudofragments.py /app/
RUN pip install --no-cache-dir --upgrade pip setuptools wheel \
 && pip install --no-cache-dir muon numpy scipy h5py anndata
CMD ["python", "h5mu_to_pseudofragments.py"]
