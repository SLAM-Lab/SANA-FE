FROM node:18-alpine AS make-sanafe

RUN apk add bash gcc g++ musl-dev python3 py3-pip python3-dev
RUN apk add make git

RUN mkdir /home/build_tutorial
RUN git clone https://github.com/SLAM-Lab/SANA-FE /home/build_tutorial/sana-fe
RUN cd /home/build_tutorial/sana-fe && make

FROM node:18-alpine AS sanafe
RUN apk add bash python3 py3-pip
RUN mkdir /tutorial
RUN python3 -m venv /tutorial/venv && . /tutorial/venv/bin/activate && pip install --upgrade pip && pip install pyyaml numpy

COPY dvs_challenge.npz /tutorial/dvs_challenge.npz
COPY --from=make-sanafe /home/build_tutorial/sana-fe/setup_tutorial.sh /tutorial/
COPY --from=make-sanafe /home/build_tutorial/sana-fe/sim /tutorial/sim
COPY --from=make-sanafe /home/build_tutorial/sana-fe/sim.py /tutorial/sim.py
COPY --from=make-sanafe /home/build_tutorial/sana-fe/arch/*.yaml /tutorial/
COPY --from=make-sanafe /home/build_tutorial/sana-fe/snn/*.net /tutorial/
COPY --from=make-sanafe /home/build_tutorial/sana-fe/scripts/challenge.py /tutorial/
RUN chmod 755 /tutorial/setup_tutorial.sh
ENTRYPOINT exec /bin/bash -c "/tutorial/setup_tutorial.sh & trap : TERM INT; sleep infinity & wait"
