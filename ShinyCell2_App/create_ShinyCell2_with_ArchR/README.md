# 📊 ShinyCell2 Docker

This project provides a **Dockerized Shiny Server** with [ShinyCell2](https://github.com/the-ouyang-lab/ShinyCell2) preinstalled.  
It allows researchers to deploy interactive single-cell visualization dashboards directly from **Seurat** or **AnnData (`.h5ad`)** objects.

---

## 🔧 Features

- Dockerized **Shiny Server** based on the official `rocker/shiny:4.5.1` image
- Preinstalled **ShinyCell2** and **ArchR**
- R packages: `Seurat`, `ggplot2`, `ggpubr`, `ggdendro`, `bslib`
- Python support: `anndata`, `scanpy` for `.h5ad` input
- Runs Shiny apps under `/srv/shiny-server/`
- Accessible from browser via `http://localhost:3838`
- Example datasets available for testing and demos

> ⚠️ **Note:** For Seurat objects, R 4.4.0 may be preferred for faster build times (~857s)

---

## 🐳 Docker Image

This Dockfile now uses a locally built base image tagged `shinycell2-local` instead of pulling from Docker Hub.
You can still reference the original image from Docker Hub if you prefer.

---

## 🚀 Usage

### 1. Build the local base image

```bash
docker build -t shinycell2-local -f ShinyCell2_App/docker/shinycell2/Dockerfile ShinyCell2_App/docker/shinycell2
```

### 2. Build this image (ArchR + RStudio)

```bash
docker build -t shinycell2-archr -f ShinyCell2_App/create_ShinyCell2_with_ArchR/Dockfile ShinyCell2_App/create_ShinyCell2_with_ArchR
```
