# mtSwirl WDL pipeline

[![DOI](https://zenodo.org/badge/448948131.svg)](https://zenodo.org/badge/latestdoi/448948131)

mtSwirl is a reference-aware pipeline for analyzing mitochondrial DNA (mtDNA)
copy number and heteroplasmy from whole-genome sequencing (WGS) data. The WDL
workflow sources are in [`inputs/wdl/`](inputs/wdl/); helper scripts and
reference resources are under [`resources/`](resources/).

This repository currently provides:

- **Single-sample v2.5:** Cromwell workflow for processing one WGS CRAM/BAM
  and producing mitochondrial calls and QC outputs.
- **Multi-sample v2.6:** serial and parallel WDL workflows for batched
  multi-sample processing. These workflows are provided as WDL sources; this
  guide's tested local run covers only the single-sample workflow.

The work in progress is validating the original single-sample workflow and
organizing reusable resources for local or HPC execution. This setup uses
four task images; it is not yet the future single-container mitoSwirl rewrite.

## Repository layout

```text
inputs/
├── config/cromwell-local-single.conf
├── json/                        # Workflow input templates
├── wdl/                         # Single- and multi-sample WDL sources
└── wgs/                          # Local sample inputs (ignored by Git)
resources/
├── cromwell/                     # Cromwell and Womtool JARs
├── docker/                       # Saved Docker image archives
├── files/                        # chrM and NUMT resources
├── grch38/                       # GRCh38 reference files
└── scripts/                      # Pipeline helper and download scripts
tests/
├── cromwell-executions/          # Execution workspaces and intermediates
├── cromwell-outputs/             # Curated final outputs
└── cromwell-workflow-logs/       # Cromwell workflow logs
```

## Run the single-sample test

The supplied offline inputs describe the NA18484 test sample. Paths in the
JSON are repository-root relative, so run validation and Cromwell from the
repository root.

### Requirements

Prepare or download inputs with internet access; the workflow itself runs
Docker tasks with networking disabled.

- Docker Engine and permission to pull and run images
- Java 17
- `curl`, `gzip`, and `md5sum`
- Google Cloud CLI (`gcloud`) to download Broad public reference files
- Substantial disk space for WGS data, references, image archives, and task
  intermediates

On Ubuntu or Debian, install Java and Docker with:

```bash
sudo apt-get update
sudo apt-get install -y openjdk-17-jre docker.io curl gzip coreutils
sudo systemctl enable --now docker
```

Install the [Google Cloud CLI](https://cloud.google.com/sdk/docs/install) if
you need to download the GRCh38 reference. The mitochondrial reference files
and supporting resources under `resources/files/` are distributed with this
package and do not need to be downloaded separately. Cromwell 92 and Womtool
92 are described below.

### Prepare inputs

The offline test JSON is `inputs/json/prepopulated_inputs_offline_single.json`.
It points to the local CRAM in `inputs/wgs/NA18484/`, GRCh38 and mitochondrial
references in `resources/`, and helper scripts in `resources/scripts/`.

Download the NA18484 CRAM and index from 1000 Genomes:

```bash
mkdir -p inputs/wgs/NA18484
cd inputs/wgs/NA18484

curl -fL --retry 5 -o NA18484.final.cram \
  'https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989343/NA18484.final.cram'
curl -fL --retry 5 -o NA18484.final.cram.crai \
  'https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989343/NA18484.final.cram.crai'

printf '%s  %s\n' \
  'e3699af246f2d14499f7ac8bd96d7ab4' NA18484.final.cram \
  'd4186393767fb0d251da9e23d07c35cc' NA18484.final.cram.crai \
  | md5sum --check
cd -
```

Place the GRCh38 FASTA, FAI, and dictionary in `resources/grch38/`. If you
need to download them, the Broad helper fetches the files from the public
reference bucket:

```bash
mkdir -p resources/grch38
cd resources/grch38
bash ../scripts/get_hg38_from_gcp.sh
cd -
```

The mitochondrial FASTA files, indices, interval lists, blacklist, and
liftover resources are already provided under `resources/files/` when using
the packaged distribution. No separate mitochondrial-reference download is
needed.

### Prepare tools

The four required task images are:

- `docker.io/rahulg603/ucsc_docker:latest`
- `docker.io/rahulg603/genomes_cloud_docker:latest`
- `docker.io/rahulg603/haplochecker_docker:latest`
- `docker.io/rahulg603/gatk_samtools_docker:latest`

Pull and save them under `resources/docker/` while online:

```bash
mkdir -p resources/docker
images=(
  docker.io/rahulg603/ucsc_docker:latest
  docker.io/rahulg603/genomes_cloud_docker:latest
  docker.io/rahulg603/haplochecker_docker:latest
  docker.io/rahulg603/gatk_samtools_docker:latest
)

for image in "${images[@]}"; do
  name="${image##*/}"
  name="${name%%:*}"
  docker pull "$image"
  docker save "$image" | gzip -c > "resources/docker/${name}.tar.gz"
done
```

Load the archives before an offline run:

```bash
for image in resources/docker/*.tar.gz; do
  docker load --input "$image"
done
```

Download Cromwell and Womtool version 92 to `resources/cromwell/`:

```bash
mkdir -p resources/cromwell
curl -fL \
  https://github.com/broadinstitute/cromwell/releases/download/92/cromwell-92.jar \
  -o resources/cromwell/cromwell-92.jar
curl -fL \
  https://github.com/broadinstitute/cromwell/releases/download/92/womtool-92.jar \
  -o resources/cromwell/womtool-92.jar
```

### Validate and run

From the repository root, load the images, validate the workflow and input
JSON, then run Cromwell:

```bash
for image in resources/docker/*.tar.gz; do
  docker load --input "$image"
done

java -jar resources/cromwell/womtool-92.jar validate \
  -i inputs/json/prepopulated_inputs_offline_single.json \
  -l inputs/wdl/fullMitoPipeline_v2_5_Single.wdl

bash resources/scripts/run_single_sample.sh
```

The config at `inputs/config/cromwell-local-single.conf` limits execution to
one task at a time, disables task networking, and sets the container hostname
to `localhost` for Hail/Spark. It routes execution directories and workflow
logs into `tests/`; the run helper copies curated outputs there after success.
The run helper waits for Cromwell to finish successfully, reads the workflow
ID from Cromwell's metadata, and automatically copies the final VCF and
diagnostic TSV into the curated output directory.

This local Cromwell backend does not enforce cloud-style WDL resource requests,
so a successful local test validates workflow execution but does not establish
HPC scheduler resource requirements.

## Outputs and cleanup

Cromwell assigns a unique workflow ID to each run. Execution data and
intermediates are retained beneath
`tests/cromwell-executions/MitochondriaPipeline/<workflow-id>/`.

The primary outputs are the final mitochondrial VCF and diagnostic statistics
TSV. The run helper copies them to
`tests/cromwell-outputs/MitochondriaPipeline/<workflow-id>/NA18484/`, leaves
the originals in the execution tree, and refuses to overwrite an existing
destination. The copy can also be run manually with
`resources/scripts/copy_single_sample_outputs.sh <workflow-id> <sample-name>`.
Workflow logs are under `tests/cromwell-workflow-logs/`.

Execution directories can be large. Only remove a run's execution directory
after confirming the workflow succeeded and keeping any intermediates or
outputs you need. Curated output copies are separate and are not removed with
the execution directory.

## Project checklist

### Completed

- [x] Run the original mtSwirl v2.5 single-sample workflow end to end on
  NA18484 with Docker task networking disabled.
- [x] Keep the test sample, references, Docker archives, Cromwell, and Womtool
  organized in the repository's `inputs/`, `resources/`, and `tests/`
  directories.
- [x] Configure the local Cromwell run and validate the single-sample WDL and
  input JSON with Womtool.
- [x] Add a run helper that automatically copies the final VCF and diagnostic
  TSV to a sample-specific output directory after a successful run.

### Remaining goals

- [ ] Build the planned single-container mitoSwirl package with its
  dependencies available locally for offline execution.
- [ ] Phase out the GATK Docker image.
- [ ] Phase out the samtools Docker image.
- [ ] Phase out the HaplotypeChecker Docker image.
- [ ] Phase out the UCSC Docker image.
- [ ] Port homoplasmy checks from R to Python.
- [ ] Port variant-bounds checks from R to Python.
- [ ] Port consensus generation from R to Python.
- [ ] Compare mitoSwirl results against the validated original v2.5 workflow
  to check output compatibility.
- [ ] Test execution on the target HPC environment and measure actual CPU,
  memory, disk, and runtime requirements; the local Cromwell run does not
  validate scheduler resource requests.
- [ ] Develop and validate the separate cross-sample merging workflow.

## Additional project context

The original mtSwirl project developed reference-aware quantification of
mtDNA copy number and heteroplasmy using WGS. The repository includes
single-sample and multi-sample WDL sources; the test instructions above focus
on validating and running the single-sample workflow.

For work using the pipeline, cite the associated paper:

> Gupta, R., Kanai, M., Durham, T.J. et al. Nuclear genetic control of mtDNA
> copy number and heteroplasmy in humans. *Nature* (2023).
> https://doi.org/10.1038/s41586-023-06426-5
