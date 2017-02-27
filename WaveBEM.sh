# Pull docker image
docker pull mathlab/wavebem:8.5.pre.4
docker run -t --rm -P -v `pwd`:/home/dealii/run:rw mathlab/wavebem:8.5.pre.4 \
    /bin/sh -c "export PATH=/home/dealii/WaveBEM-travis-to-latest/build/:\$PATH; cd /home/dealii/run; WaveBEM  $@"
