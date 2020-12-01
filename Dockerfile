# Using OpenJDK 8
FROM broadinstitute/gatk:latest

RUN apt-get update \
    && apt-get install -y git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ADD . /discvr-build

RUN cd /discvr-build \
    && ./gradlew assemble \
    && ./gradlew installDist \
    && ./gradlew shadowJar \
    && ls build/libs/ \
    && mv build/libs/DISCVRSeq-*.jar /DISCVRSeq.jar \
    && cd / \
    && rm -Rf /discvr-build

ENTRYPOINT ["java", "-jar", "DISCVRSeq.jar"]
