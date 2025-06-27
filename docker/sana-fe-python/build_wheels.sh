#!/bin/bash

ROOTLESS_OR_NOT=$(docker info | grep -i "rootless")
if [[ -z "$ROOTLESS_OR_NOT" ]]; then
    SOCKET_PATH="/var/run/docker.sock"
    echo "You are not running Docker in rootless mode. Using default socket path: $SOCKET_PATH"
else
    SOCKET_PATH="$XDG_RUNTIME_DIR/docker.sock"
    echo "You are running Docker in rootless mode. Using socket path: $SOCKET_PATH"
fi

if [[ ! -S "$SOCKET_PATH" ]]; then
    echo "Docker socket not found at $SOCKET_PATH. Make sure that Docker is running and accessible."
    exit 1
else
    echo "Using Docker socket at $SOCKET_PATH."
fi

docker build -t sana-fe-python .

if [[ $? -ne 0 ]]; then
    echo "Docker build failed."
    exit 1
fi

mkdir -p wheelhouse

docker run --rm -v "$(pwd)/wheelhouse:/app/wheelhouse" -v "$SOCKET_PATH:/var/run/docker.sock" sana-fe-python