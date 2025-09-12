# microSplit_scRNAseq

## Setup

Run `setup.sh` to set default paths

```
source setup.sh
```

To obtain container locally (if not already downloaded) (Note for covert lab: we already have it on OAK)
```
singularity pull docker://scr.svc.stanford.edu/khoang99/containers/scrna:latest
# set RUN environment variable to execute singularity run
export RUN="singularity run scrna_latest.sif"
```

## Run Jupyter Notebook

To run Jupyter Notebook with container, first run jupyter server on your compute node (Sherlock), ideally in a new `screen` (e.g., `screen -S notebook`)

```
# port can be set to any number from 20000 to 30000, pick a random one to not overlap with other users
$RUN jupyter notebook --port=23232 --no-browser --ip=127.0.0.1
```

Now, on your local computer, login to Sherlock with ssh forwarding

```
# modify port number according to what you chose above
ssh [USER]@login.sherlock.stanford.edu -L 23232:localhost:23232 
# When inside login node, ssh port forward to compute node (modify sh02-02n26 according to compute node you start the jupyter server)
ssh -N -L 23232:localhost:23232 sh02-02n26
```

On local machine, paste http://localhost:23232/ (adjust port accordingly) into your browser, enter Sherlock password to login Jupyter notebook
