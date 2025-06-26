FROM ubuntu:latest AS build

WORKDIR /opt/app
COPY . .

RUN apt-get update && apt-get install -y python3 python3-pip python3-dev make gcc g++ \
    git flex bison cmake pybind11-dev && \
    apt-get clean

RUN mkdir build && \
    cd build && \
    cmake -DPYTHON_BUILD_ENABLED=OFF -DSTANDALONE_BUILD_ENABLED=ON .. && \
    make

ENTRYPOINT [ "/opt/app/build/sim" ]