# ATAC2GRN Project README

This is the Github for ATAC2GRN, a validated pipeline for running chromatin profiling assays. Here you'll find code that was used to produce the figures in the paper as well as pipelines and tools for validating your own pipelines.

## How to Get Started

If you're looking to take bits and pieces of this code for your own use, skip forward to that folder's discussion. Otherwise, here are two options to set up your machine with the appropriate prerequisites.

### Apptainer File

Now with a new Apptainer image and definition file for immediate rollout on your high-performance computing cluster of choice! Download the new .sif file [here](https://drive.google.com/file/d/1UCvvgN3Bs85otx_AXmIazv1GznzPgU4C/view?usp=sharing).

For more information on how to use Apptainer in your system, visit [here](https://apptainer.org/docs/user/latest/quick_start.html).

### Using VirtualBox and a VDI File

If you'd like to build a machine in VirtualBox to run these pipelines, download the .vdi file from [our Google Drive](https://drive.google.com/open?id=1j-sO0CjyK-u95Y2ZPQIDz5zqH2KnvAnd).

Unpack the file with

```
tar -xvzf ATAC2GRN_VDI.tar.gz
```

Inside the folder will be a .vdi file, and this can be loaded into VirtualBox directly as an Ubuntu 16.04 system with all of the dependencies for this pipelines installed. We'd recommend allocating sizable hard drive space. We used 400 GB, though as a rule of thumb you should have at least 250% as much space as it takes to accomodate your raw data. Additional make sure you have enough RAM (>= 16 GB) so that running the pipeline doesn't take multiple days. A CPU that can run at least four threads will work with the code.

## The Folders

There are three folders in this project: the code used to generate figures; the code used for pipelines in both bash and Snakemake; and the tools used to assess pipeline recapitulation of ChIP-seq.

### Code to Generate Figures

This code holds unit-tested code for figure generation in the original paper.

### Pipes

This folder contains all the pipeline code used in the paper, including the bash code used to test all 4560 pipelines as well as the Snakemake alternative for final pipelines. Interested users are provided with both, to use depending on their personal preference.

To make best use of these tools, please open the Pipes folder and follow the README.md instructions on how to use these finalized pipelines. 

### ROCtool

This folder contains the code necessary to test a set of footprints against known ChIP-seq data. This ChIP-seq data is specific to GM12878; a new directory of ChIP-seq data needs to be made if a user produces footprint data from a different pipeline.

Information on how to use this tool is provided in the ROCtool folder. The tool uses the .bed files of all footprints produced by a finalized pipeline.