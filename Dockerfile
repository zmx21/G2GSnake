FROM rocker/r-ver:4.1.0

RUN R -e "install.packages('phylobase',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('adephylo',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('filematrix',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('seqinr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('purr',dependencies=TRUE, repos='http://cran.rstudio.com/')"

WORKDIR /home/

RUN mkdir ./results
RUN mkdir ../raw_data/
RUN mkdir ./scripts

COPY plink2 ./plink2

COPY gcta-1.94.1-linux-kernel-3-x86_64.zip gcta-1.94.1-linux-kernel-3-x86_64.zip
RUN unzip gcta-1.94.1-linux-kernel-3-x86_64
RUN cp gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 ./gcta64

COPY bcftools-1.17.tar.bz2 bcftools-1.17.tar.bz2

RUN tar -xf bcftools-1.17.tar.bz2
RUN apt update
RUN apt install zlib1g zlib1g-dev libbz2-dev liblzma-dev
WORKDIR ./bcftools-1.17
RUN ./configure
RUN make 
RUN make install 

WORKDIR /home/
