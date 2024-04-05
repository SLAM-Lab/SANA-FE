FROM node:18-alpine

RUN apk add bash gcc g++ musl-dev python3 py3-pip python3-dev
RUN apk add make git

RUN mkdir /home/tutorial
RUN python3 -m venv /home/tutorial/venv && . /home/tutorial/venv/bin/activate && pip install --upgrade pip && pip install pyyaml numpy
COPY dvs_challenge.npz /home/tutorial/dvs_challenge.npz

