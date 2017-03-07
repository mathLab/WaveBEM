FROM mathlab/deal2lkit:8.5.pre.4-debugrelease

MAINTAINER luca.heltai@gmail.com
ARG BUILD_TYPE=DebugRelease

ARG VER=master

RUN \
    git clone https://github.com/mathLab/WaveBEM/ WaveBEM-$VER && \
    cd WaveBEM-$VER && git checkout $VER && \
    mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
	  -DENABLE_DOCUMENTATION=OFF \
	  -GNinja \
          ../ && \
    ninja -j4 

WORKDIR $HOME/WaveBEM-$VER/build
