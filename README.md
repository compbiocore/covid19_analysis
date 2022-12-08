# Covid19 Docker Container

Template for analyses repositories. For more information see https://compbiocore-brown.slab.com/posts/data-organisation-for-analysis-repos-fdi2cddd. Folders that should be present in all such repositories are:

 * **metadata:** contains ```Dockerfile``` and a fasta file for testing ```ri_sequences.fa``` comprised of 50 covid genome sequences downloaded from NCBI
 * **scripts:** contains ```run.sh``` 

Build the docker container using:

```
docker build --platform linux/x86_64 -t covid19 .
```

To initiate and instance of the container run:
```
docker run -it covid19:latest bash
```

# CBC Project Information

```
title: Covid19 docker container
tags:
analysts:
git_repo_url:
resources_used: Pangolin, Nextclade, Nextalign
summary: 
project_id:
```
