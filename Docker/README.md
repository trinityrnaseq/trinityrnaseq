
# build Trinity docker

    docker build -t trinityrnaseq/trinityrnaseq:${RELEASE_TAG} .


# run Trinity via docker

    docker run trinityrnaseq/trinityrnaseq:${RELEASE_TAG} Trinity
