# Install R version 4.0.5
FROM wfma888/plaqviewmaster:latestDec10

# Install Shiny server
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt) && \
    wget "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb
    

# Copy configuration files into the Docker image
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
RUN rm -rf /srv/shiny-server/*
# COPY /* /srv/shiny-server/

# Get the app code
RUN git clone https://github.com/MillerLab-CPHG/PlaqView.git
RUN cp -r PlaqView/* /srv/shiny-server/

# Make the ShinyApp available at port 80
EXPOSE 80
WORKDIR /srv/shiny-server
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');shiny::runApp('app.R')"

#RUN chown shiny.shiny /usr/bin/shiny-server.sh && chmod 755 /usr/bin/shiny-server.sh

# Run the server setup script
#CMD ["/usr/bin/shiny-server.sh"]