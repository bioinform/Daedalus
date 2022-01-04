FROM continuumio/miniconda3:4.10.3

RUN apt-get update \
    && apt-get install -y build-essential procps

#RUN useradd -ms /bin/bash  daedalus
#USER daedalus
#WORKDIR /home/daedalus

COPY environment.yml .
COPY install.sh .

RUN conda env create -f environment.yml

RUN chmod 755 /root
RUN echo "conda activate Daedalus" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

RUN export GIT_SSL_NO_VERIFY=1 && bash install.sh
