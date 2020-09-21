FROM trinityrnaseq/trinityrnaseq:latest
SHELL ["/bin/bash", "-c"]

ADD ./parse_commands.py /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

RUN apt-get update && apt-get install --no-install-recommends -y parallel \
    && chmod a+rx /software/monitor_script.sh


