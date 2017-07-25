FROM python:2.7

WORKDIR /usr/src/app

RUN apt-get update \
    && apt-get install -y vim libopenmpi-dev openmpi-bin gfortran
COPY requirements.txt ./

## Need to get numpy prior to all requirements...
RUN pip install --upgrade pip setuptools wheel numpy

## Install requirements
RUN pip install -r requirements.txt

COPY . .

ENV PYTHONPATH $PYTHONPATH:$PWD

RUN make

CMD [ "./coverage_and_test.sh"]
