# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-2.2.0

RUN apt-get update \
    && apt-get install -y git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ADD . /discvr-build

RUN cd /discvr-build \
    && export GIT_TAG=$(git describe --tags) \
    && echo $GIT_TAG \
    && ./gradlew assemble \
    && ./gradlew installDist \
    && ./gradlew shadowJar -Drelease=true \
    && ls build/libs/ \
    && mv build/libs/DISCVRSeq-${GIT_TAG}.jar /DISCVRSeq.jar \
    && cd / \
    && rm -Rf /discvr-build

