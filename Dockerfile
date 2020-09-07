# Source Image
FROM centos:7

RUN yum update -y && yum install -y git gcc-c++ zlib-devel make gzip boost-devel boost-static

COPY ./ /duohmm/
WORKDIR /duohmm/
RUN BOOST_ROOT=/usr/ make -j4 
RUN make test
