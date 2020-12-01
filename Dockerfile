# Using OpenJDK 8
FROM broadinstitute/gatk:latest

RUN apt-get update \
    && apt-get install -y git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ADD . /discvr-build

RUN cd /discvr-build \
    && export GIT_TAG=$(git describe --tags --abbrev=0) \
    && echo $GIT_TAG \
    && ./gradlew assemble \
    && ./gradlew installDist \
    && ./gradlew shadowJar -Drelease=true \
    && ls build/libs/ \
    && mv build/libs/DISCVRSeq-${GIT_TAG}.jar /DISCVRSeq.jar \
    && cd / \
    && rm -Rf /discvr-build

ENTRYPOINT ["java", "-jar", "DISCVRSeq.jar"]
