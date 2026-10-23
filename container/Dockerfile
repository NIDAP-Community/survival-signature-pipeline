# syntax=docker/dockerfile:1
#
# Portable OCI image for the maintained survival-analysis CLI. The image is
# intentionally independent of any AWS, database, or OncoChat credentials:
# callers mount analysis inputs at /data and receive run artifacts at /output.

ARG R_VERSION=4.5.1
FROM rocker/r-ver:${R_VERSION}

ARG VCS_REF=unknown
ARG BUILD_DATE=unknown
LABEL org.opencontainers.image.title="survival-signature-pipeline" \
      org.opencontainers.image.description="Reusable R survival analysis pipeline" \
      org.opencontainers.image.source="https://github.com/NIDAP-Community/survival-signature-pipeline" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}"

WORKDIR /app

# These packages cover the compiled and plotting dependencies pinned in
# renv.lock, including curl/xml, text shaping, raster plotting, and Rcpp-based
# survival-model packages. Keep this list explicit so image builds are
# reproducible and diagnose missing system libraries clearly.
RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        g++ \
        gfortran \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgit2-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff-dev \
        libxml2-dev \
        make \
    && rm -rf /var/lib/apt/lists/*

# Restore dependencies before copying the application source. This keeps the
# expensive package-install layer cacheable when only R code or documentation
# changes. --vanilla avoids the repository's local interactive activation path.
COPY renv.lock ./
COPY renv/activate.R renv/settings.json ./renv/
RUN Rscript --vanilla -e 'install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore(project = "/app", prompt = FALSE)'

COPY . ./
COPY container/entrypoint.sh /usr/local/bin/survival-pipeline
RUN chmod 0755 /usr/local/bin/survival-pipeline \
    && mkdir --parents /data /output \
    && useradd --create-home --uid 10001 --shell /usr/sbin/nologin survival \
    && chown --recursive survival:survival /app /data /output

# Do not evaluate user profiles in a noninteractive run. The entrypoint
# explicitly resolves the restored renv library and passes it to every child
# Rscript process launched by entrypoint.R.
ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    R_PROFILE_USER=/dev/null

USER survival
ENTRYPOINT ["/usr/local/bin/survival-pipeline"]
CMD ["list"]
