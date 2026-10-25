# WDL workflows and inputs

Several reference file and Docker image inputs are required for proper and reproducible functioning of the mtSwirl pipeline.

## Code

Supporting scripts and reference resources are stored in `resources/scripts`
and `resources/files`. The input templates are in `inputs/json`; each helper
file should be supplied to the workflow input with the corresponding field.

## Docker images

The Docker images correspond to the following workflow inputs:

- `ucsc_docker`: rahulg603/ucsc_docker
- `genomes_cloud_docker`: rahulg603/genomes_cloud_docker
- `haplochecker_docker`: rahulg603/haplochecker_docker
- `gatk_samtools_docker`: rahulg603/gatk_samtools_docker

**Please note that Docker Hub places pull limits making these links unusable for large scale parallelized workflows.** We recommend using `docker pull` locally to obtain the image, then pushing the image to either GCP artifact registry (see [here](https://cloud.google.com/artifact-registry/docs/docker/store-docker-container-images) for a quickstart guide) for Terra/Cromwell-based executions or to DNANexus (see [here](https://documentation.dnanexus.com/developer/apps/dependency-management/using-docker-images)).

## Multi-sample execution on dx

There are some additional notes required for execution on DNANexus. Once you have ensured that you have proper access to the WGS files of interest on the platform, dxCompiler must first be used to upload the workflow:

```
bash resources/scripts/run-dxcompiler-docker-multi compile inputs/wdl/fullMitoPipeline_v2_6_Multi.wdl -project PROJECT_ID -folder DESIRED_DX_PATH -f
```

**UPDATE:** as of v2 of this pipeline, we have developed a new version which uses mutli-core processing per-VM to greatly boost compute speed for time-intensive tasks. To run the parallelized pipeline, use:

```
bash resources/scripts/run-dxcompiler-docker-multi compile inputs/wdl/parallel_fullMitoPipeline_v2_6_Multi.wdl -project PROJECT_ID -folder DESIRED_DX_PATH -f
```

We have validated that this version of the pipeline produces identical results to the serial version. To take advantage of this pipeline, you must use larger VMs. You can see some example sizes in `inputs/json/prepopulated_inputs_dx_parallel.json`,
which contains the new variables that need to be used with this WDL.

For DNANexus, we have several recommendations when building input schemas for batch job submission, based on `inputs/json/prepopulated_inputs_dx.json` or `inputs/json/prepopulated_inputs_dx_parallel.json`:

1. Please do not use docker.io paths for large scale workflows, and instead upload images to your dx project (see [here](https://documentation.dnanexus.com/developer/apps/dependency-management/using-docker-images)).
2. We recommend saving all GCP-hosted files (as listed in `inputs/json/prepopulated_inputs_dx.json`) into your own project directory on the dx platform, and changing paths to the corresponding [DNANexus Link](https://documentation.dnanexus.com/developer/api/running-analyses/job-input-and-output).
3. To run this pipeline, three additional items must be added to the `prepopulated_inputs_dx.json` input used for each workflow run:
   1. `stage-common.sample_name`: an array of sample names
   2. `stage-common.wgs_aligned_input_bam_or_cram`: an array of DNANexus links corresponding to WGS CRAM files
   3. `stage-common.wgs_aligned_input_bam_or_cram_index`: an array of DNANexus links corresponding to WGS CRAI files

To add these additional fields, we recommend batching the list of all CRAMs/CRAIs into 18 samples (fullMitoPipeline) or 40 samples (fullMitoSamples). Each batch is processed on one machine at a time either serially or in parallel respectively.