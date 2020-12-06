# Using OpenJDK 8
FROM broadinstitute/gatk:latest

# See: https://stackoverflow.com/questions/44331836/apt-get-install-tzdata-noninteractive
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y wget ncbi-blast+ build-essential g++ cmake git-all \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

#Primer3
RUN cd / \
    && wget https://github.com/primer3-org/primer3/archive/v2.5.0.tar.gz \
    && tar -xf v2.5.0.tar.gz \
    && rm v2.5.0.tar.gz \
    && cd /primer3-2.5.0/src \
    && make

ENV PATH="/primer3-2.5.0/src:${PATH}"

ADD . /discvr-build

RUN cd /discvr-build \
    && ./gradlew assemble \
    && ./gradlew shadowJar \
    && ./gradlew copyShadowJar \
    && mv build/jars/DISCVRSeq-*.jar /DISCVRSeq.jar \
    && cd / \
    && rm -Rf /discvr-build

# Santity check:
#RUN which blastn
#RUN which primer3_core
#stat /DISCVRSeq.jar

ENTRYPOINT ["java", "-jar", "/DISCVRSeq.jar"]
