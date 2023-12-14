FROM rocker/verse:4.2.2

COPY . frasyr

RUN R -e "devtools::install_local('frasyr')"

ENV EDITOR_FOCUS_DIR="/home/rstudio/frasyr"

RUN mkdir -p $EDITOR_FOCUS_DIR

USER 0
RUN echo "session-default-working-dir=$EDITOR_FOCUS_DIR" >> /etc/rstudio/rsession.conf && \
    echo "session-default-new-project-dir=$EDITOR_FOCUS_DIR" >> /etc/rstudio/rsession.conf
