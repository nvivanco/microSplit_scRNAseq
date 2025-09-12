# microSplit_scRNAseq

## Setup

Run `setup.sh` to set default paths

```
source setup.sh
```

To obtain container locally
```
singularity pull docker://scr.svc.stanford.edu/khoang99/containers/scrna:latest
export RUN="singularity run scrna_latest.sif"
```
