#!/bin/bash

IMAGE_NAME="sana-fe"

docker build -t $IMAGE_NAME .
if [ $? -eq 0 ]; then
    echo "Docker image '$IMAGE_NAME' built successfully."
else
    echo "Failed to build Docker image '$IMAGE_NAME'."
    exit 1
fi

docker run --rm -v "$(pwd)/arch":/data/arch -v "$(pwd)/snn":/data/snn $IMAGE_NAME "/data/$1" "/data/$2" "${@:3}"