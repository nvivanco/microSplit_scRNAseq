# microSplit_scRNAseq

## Setup

Run `setup.sh` to set default paths

```
source setup.sh
```

To obtain container locally (if not already downloaded)
```
singularity pull docker://scr.svc.stanford.edu/khoang99/containers/scrna:latest
# set RUN environment variable to execute singularity run
export RUN="singularity run scrna_latest.sif"
```
