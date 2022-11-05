FROM ubuntu:bionic

RUN apt-get -y update && apt-get -y upgrade
# The following is necessary to avoid an interactive prompt when installing r-base
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
# instructions here: https://www.rstudio.com/products/shiny/download-server/ubuntu/

#
# from https://cloud.r-project.org/bin/linux/ubuntu/
#
# update indices
apt update -qq
# install two helper packages we need
apt install -y --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.2 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran42/"



RUN apt-get install -y r-base r-base-dev
RUN apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev

RUN Rscript -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN apt-get install -y gdebi-core wget
RUN wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.19.995-amd64.deb
RUN gdebi --n shiny-server-1.5.19.995-amd64.deb

# remove the default landing page and link to sample app's
RUN rm /srv/shiny-server/*

# overwrite the default config with our modified copy
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf 
RUN chmod 777 /etc/shiny-server/shiny-server.conf

# This is the app folder specified in shiny-server.conf
RUN mkdir -p /srv/shiny-server/app

# make the installation folder and library folder accessible to the 'shiny' user
RUN chmod -R 777 /srv/shiny-server/
RUN chmod -R 777 /usr/local/lib/R/site-library
RUN chmod -R 777 /var/lib/shiny-server

# This is where the app' will be installed
WORKDIR /srv/shiny-server/app

# Run the server as the 'shiny' user
USER shiny

# Send application logs to stderr
ENV SHINY_LOG_STDERR=1
ENV SHINY_LOG_LEVEL=TRACE


#
# Above is generic Shiny server set up
# Below is specific to this app'
#

#FROM sagebionetworks/shiny-base:release-1.0
USER root
RUN apt-get install -y default-jdk
USER shiny
# This is the expected application installation folder
WORKDIR /srv/shiny-server/app
COPY --chown=shiny ./ ./
# renv restore
RUN Rscript -e "install.packages(c('renv'), repos='http://cran.rstudio.com/'); renv::restore()"
RUN Rscript -e "install.packages(c('conflicted','rjson'), repos='http://cran.rstudio.com/')"
