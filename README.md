[![Build Status](https://api.travis-ci.com/BimberLab/DISCVRSeq.svg)](https://travis-ci.com/BimberLab/DISCVRSeq)

## Overview
DISCVR-seq Toolkit is a diverse collection of tools for working with sequencing data, developed and maintained by the Bimber Lab, built using the GATK4 engine. The set of tools is analogous to GATK or Picard.  A description of all software produced by the Bimber Lab can be found [here](https://bimberlab.github.io).  

## Documentation
[Please view our documentation](https://bimberlab.github.io/DISCVRSeq/) for more information about the set of tools and usage.

## Getting Started
DISCVR-seq Toolkit is a java program distributed as a single JAR.  You can download the latest JAR from our [release page](https://github.com/BimberLab/DISCVRSeq/releases).  Running tools is analogous to GATK4.  

While we recommend [our documentation](https://bimberlab.github.io/DISCVRSeq/) to learn about available tools and options, one can also view a list of tools and/or arguments from the command line:

```

# List available tools:
java -jar DISCVRseq.jar --list 

# View arguments for a specific tool (VariantQC in this example):
java -jar DISCVRseq.jar VariantQC --help

```

## Docker

By popular demand, DISCVR-seq releases  [will now be pushed to Docker Hub](https://hub.docker.com/repository/docker/bbimber/discvrseq).  We recommend using a specific release, which you can do using tags:

docker pull bbimber/discvrseq:release1.11