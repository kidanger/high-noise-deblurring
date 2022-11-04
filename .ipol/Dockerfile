FROM python:3.5

RUN apt-get update && apt-get install -y build-essential libpng-dev libjpeg-dev libtiff-dev libeigen3-dev cmake libfftw3-dev bc && rm -fr /var/cache/apt

# copy the requirements.txt and install python packages
COPY ffdnet-denoiser/requirements.txt requirements.txt
#RUN pip3 install --no-cache-dir torch==0.4.0 -f https://download.pytorch.org/whl/cpu/torch_stable.html
RUN pip3 install --no-cache-dir -r requirements.txt && rm requirements.txt

RUN wget -q https://www.ipol.im/pub/art/2018/211/imgutils-master.zip && unzip imgutils-master.zip && cd imgutils-master && cmake . && make -j
RUN wget -q https://www.ipol.im/pub/art/2018/211/imscript.zip && unzip imscript.zip && make -C imscript-master -j

# copy the code to $bin
ENV bin /workdir/bin/
RUN mkdir -p $bin
WORKDIR $bin
COPY . .

RUN make -j -C kernel-estimation && make -j -C ponomarenko-noise-estimation && mkdir nonblind-deconvolution/build && cd nonblind-deconvolution/build && cmake .. && make -j

# the execution will happen in the folder /workdir/exec
# it will be created by IPOL

# some QoL tweaks
ENV PYTHONDONTWRITEBYTECODE 1
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION python
ENV PATH $bin:$PATH:/imgutils-master:/imscript-master/bin

# $HOME is writable by the user `ipol`, but 
ENV HOME /home/ipol
# chmod 777 so that any user can use the HOME, in case the docker is run with -u 1001:1001
RUN groupadd -g 1000 ipol && useradd -m -u 1000 -g 1000 ipol -d $HOME && chmod -R 777 $HOME
USER ipol
