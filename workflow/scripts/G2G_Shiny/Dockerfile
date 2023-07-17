FROM rocker/shiny:4.2.1
ENV RENV_CONFIG_REPOS_OVERRIDE https://packagemanager.rstudio.com/cran/latest

RUN apt-get update -qq && \ 
  apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libicu-dev \
    libssl-dev \
    make \
    pandoc \
    zlib1g-dev && \
  apt-get clean && \ 
  rm -rf /var/lib/apt/lists/*

ENV _R_SHLIB_STRIP_=true

COPY shiny_renv.lock renv.lock
RUN Rscript -e "install.packages('renv')"
RUN Rscript -e "renv::restore()"

COPY app.R /srv/shiny-server/

RUN mkdir /home/results/

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
