FROM python:3.6

MAINTAINER Julien Peloton <j.peloton@sussex.ac.uk>

WORKDIR /usr/src/app

## Install vim, MPI, and fortran compiler
RUN apt-get update \
    && apt-get install -y vim libopenmpi-dev openmpi-bin gfortran

## Copy repo files
COPY . .

## Need to get numpy prior to all requirements (distutils)
RUN pip install --upgrade pip setuptools wheel numpy==1.14.5 \
    && pip install -r requirements.txt \
    && pip install pyslalib>=1.0.4

ENV PYTHONPATH $PYTHONPATH:$PWD

## Compile fortran sources
RUN make

CMD ["./coverage_and_test.sh"]
