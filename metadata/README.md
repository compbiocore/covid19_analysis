

Contents of `2_metadata` reflect the files in `/gpfs/data/ris3/2_metadata` which are essential for pangoling, nextclade, and nextalign analyses.

The ```Dockerfile``` is the initial file used for the creation of the container in which analyses are run.

If you want to rebuild the docker container using:

```
docker build --platform linux/x86_64 -t covid19 .
```

To initiate and instance of the container run:
```
docker run -it covid19:latest bash
```

The image can also pulled directly onto oscar as a singularity image using ```singularity pull docker://ericsalomaki/covid12162022:latest```