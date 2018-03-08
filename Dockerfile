############################################################
# Dockerfile to build GLASSgo Containers
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu

# File Author / Maintainer
MAINTAINER Steffen C. Lott

# Update the sources list
RUN apt-get update
RUN apt-get install -y gcc
RUN apt-get install -y python python-pip
RUN pip install --upgrade pip
RUN apt-get install -y python3 wget python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install numpy
RUN pip3 install biopython
RUN apt-get install -y ncbi-blast+
RUN apt-get install -y clustalo
RUN apt-get install -y bioperl
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:j-4/vienna-rna
RUN apt-get update
RUN apt-get install -y vienna-rna

# Copy the application folder inside the container
ADD /GLASSGO /bin/GLASSgo/

# Set the default directory where CMD will execute
WORKDIR /bin/GLASSgo/
