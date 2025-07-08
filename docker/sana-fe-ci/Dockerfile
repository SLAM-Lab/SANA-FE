FROM ubuntu:latest AS build

WORKDIR /opt/app
COPY . .

RUN apt-get update && apt-get install -y python3 python3-pip python3-dev make gcc g++ \
    git flex bison cmake pybind11-dev ruby clang valgrind python3-venv && \
    apt-get clean

RUN cp ci/* .

RUN export GCC_PATH=$(which g++) && \
    export CLANG_PATH=$(which clang)

CMD ["ruby", "run.rb"]