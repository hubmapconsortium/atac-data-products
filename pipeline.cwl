#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for concatenating sc-ATAC-seq datasets into MuData object

requirements:
  ScatterFeatureRequirement: {}

inputs: 
    data_directory:
        label: "Path to directory containing cell by gene and cell by bin files"
        type: Directory
    
    uuids_file:
        label: "Path to a file containing a list of uuids and other metadata for the dataset to be indexed"
        type: File
    
    tissue:
        label: "Two letter tissue type code"
        type: string?

outputs:
    mudata_file:
        outputSource: concatenate/mudata_file

steps:

    - id: concatenate
      in: 
        - id: data_directory
          source: data_directory
        - id: uuids_file
          source: uuids_file
        - id: tissue
          source: tissue
    
      out:
        - mudata_file
      run: steps/concatenate.cwl
      label: "Concatenates h5ad files in directory"
