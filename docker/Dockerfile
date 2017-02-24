FROM ubuntu:latest
MAINTAINER Christian Schafmeister <meister@temple.edu>
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y gcc g++ llvm clang cmake libgc-dev libgmp-dev binutils-gold binutils-dev zlib1g-dev libncurses-dev libboost-filesystem-dev libboost-regex-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-iostreams-dev
RUN apt-get install -y git
RUN git clone -b master https://github.com/drmeister/externals-clasp ~/externals-clasp
RUN cd ~/externals-clasp && make
RUN apt-get install -y sbcl
RUN git clone -b testing https://github.com/drmeister/clasp ~/clasp
RUN cd ~/clasp && echo "EXTERNALS_CLASP_DIR = '$HOME/externals-clasp'" > wscript.config
RUN cd ~/clasp/extensions && git clone -b dev https://github.com/drmeister/cando.git cando
RUN cd ~/clasp && ./waf update_submodules
RUN cd ~/clasp && ./waf configure
RUN cd ~/clasp && ./waf -j $(nproc) build_cboehm
RUN cd ~/clasp/build && ./clasp
RUN git clone https://github.com/slime/slime ~/slime
RUN git clone https://github.com/quicklisp/quicklisp-client.git ~/quicklisp
RUN /root/clasp/build/clasp -e "(load \"/root/slime/swank-loader.lisp\")" -e "(swank-loader:init :delete nil :reload nil :load-contribs nil)" -e "(quit)"
ENTRYPOINT /bin/bash