# syntax=docker/dockerfile:1
#
# Portable OCI image for the maintained survival-analysis CLI. The image is
# intentionally independent of any AWS, database, or OncoChat credentials:
# callers mount analysis inputs at /data and receive run artifacts at /output.

ARG R_VERSION=4.5.1
FROM rocker/r-ver:${R_VERSION}

# l2psupp supplies updategenes() for auditable normalization of user-entered
# human gene symbols. Pin the upstream repository revision instead of a mutable
# branch so the container's resolver database is reproducible.
ARG L2P_COMMIT=ab389bc978441f1893414da7a07e803bea971205

WORKDIR /app

# These packages cover verified HTTPS retrieval plus the compiled and plotting
# dependencies pinned in renv.lock, including curl/xml, text shaping, raster
# plotting, and Rcpp-based survival-model packages. Keep this list explicit so
# image builds are reproducible and diagnose missing system libraries clearly.
RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        ca-certificates \
        cmake \
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
        libx11-dev \
        libxml2-dev \
        make \
        pandoc \
    && rm -rf /var/lib/apt/lists/*

# Restore dependencies before copying the application source. This keeps the
# expensive package-install layer cacheable when only R code or documentation
# changes. --vanilla avoids the repository's local interactive activation path.
COPY renv.lock ./
COPY renv/activate.R renv/settings.json ./renv/
# A caller on an HTTPS-inspecting network may pass its trust bundle as the
# optional build_ca secret. It is used only while retrieving dependencies and
# is never copied into an image layer or repository artifact.
# Retain the renv package cache outside the image layer so locally interrupted
# builds resume completed package restores. Copy packages into the project
# library instead of linking to that ephemeral cache so the final image is
# self-contained at runtime.
RUN --mount=type=secret,id=build_ca,required=false \
    --mount=type=cache,target=/root/.cache/R/renv \
    if [ -s /run/secrets/build_ca ]; then export CURL_CA_BUNDLE=/run/secrets/build_ca; fi; \
    export RENV_CONFIG_CACHE_SYMLINKS=FALSE; \
    Rscript --vanilla -e 'install.packages("renv", repos = "https://cloud.r-project.org"); library(renv); project_library <- renv::paths$library(project = "/app"); renv::restore(project = "/app", library = project_library, prompt = FALSE)'
RUN --mount=type=secret,id=build_ca,required=false \
    if [ -s /run/secrets/build_ca ]; then export CURL_CA_BUNDLE=/run/secrets/build_ca; fi; \
    Rscript --vanilla -e "install.packages('https://github.com/CCBR/l2p/raw/${L2P_COMMIT}/l2psupp_0.0-14.tar.gz', repos = NULL); stopifnot(as.character(utils::packageVersion('l2psupp')) == '0.0.14'); stopifnot(l2psupp::updategenes('CDCA1')\$newname[[1]] == 'NUF2')"

COPY . ./
COPY container/entrypoint.sh /usr/local/bin/survival-pipeline

# Put changing provenance labels after expensive dependency layers so a new
# local build timestamp does not force a full dependency restore.
ARG VCS_REF=unknown
ARG BUILD_DATE=unknown
LABEL org.opencontainers.image.title="survival-signature-pipeline" \
      org.opencontainers.image.description="Reusable R survival analysis pipeline" \
      org.opencontainers.image.source="https://github.com/NIDAP-Community/survival-signature-pipeline" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}"

RUN chmod 0755 /usr/local/bin/survival-pipeline \
    && mkdir -p /data /output \
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
