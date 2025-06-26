FROM ubuntu:latest

# Add Docker's official GPG key (for docker-ce-cli installation):
RUN apt-get update && apt-get install -y \
    ca-certificates \
    curl \
    gnupg \
    lsb-release && \
    mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc && \
    chmod a+r /etc/apt/keyrings/docker.asc && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
    $(. /etc/os-release && echo ${UBUNTU_CODENAME:-focal}) stable" > /etc/apt/sources.list.d/docker.list && \
    apt-get update

RUN apt update && apt install -y \
    python3 \
    python3-pip \
    python3-dev \
    docker-ce-cli

WORKDIR /app
COPY . /app

RUN pip3 install --break-system-packages cibuildwheel pybind11

# workaround for cibuildwheel issue with flex
CMD ["bash", "-c", "cibuildwheel --output-dir wheelhouse"]