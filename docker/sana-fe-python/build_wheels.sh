#!/bin/bash

SOCKET_PATH="/var/run/docker.sock"

if [[ ! -S "$SOCKET_PATH" ]]; then
    echo "Docker socket not found at $SOCKET_PATH. Make sure that Docker is running and accessible."
    exit 1
fi

docker build -t sana-fe-python .

if [[ $? -ne 0 ]]; then
    echo "Docker build failed."
    exit 1
fi

mkdir -p wheelhouse

docker run --rm -v "$(pwd)/wheelhouse:/app/wheelhouse" -v "$SOCKET_PATH:$SOCKET_PATH" sana-fe-python