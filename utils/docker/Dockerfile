# syntax=docker/dockerfile:1.3

# Note: two-stage building is disabled for now because of the dual nature of
# varfish-server-worker right now as both CLI and lib; we would have to use features to
# untable the two.

# Based on https://levelup.gitconnected.com/1940db638a6c
#
# We don't do cross compilation at the moment but build the dependencies first
# anyway to get the ability to increment.

# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

# Pinning Rust version for now because of this issue:
#
# - https://github.com/rust-lang/rust/issues/95926
FROM rust:1-bookworm AS builder

# Build dependencies first.
#
# Install dependencies for compilation of C code (e.g., rocksdb).
RUN apt-get update && \
    apt-get install -y clang
# Add the needed Cargo components.
RUN rustup component add rustfmt
# Install build dependency `protoc`.
COPY utils/install-protoc.sh /tmp
RUN PREFIX=/usr/local bash /tmp/install-protoc.sh

# Now for the two-step building.
#
# Set initial workdir.
WORKDIR /usr/src
# Create blank project.
RUN USER=root cargo new varfish-server-worker
# We want dependencies cached, so copy those first.
COPY Cargo.toml Cargo.lock /usr/src/varfish-server-worker/
# Set the working directory.
WORKDIR /usr/src/varfish-server-worker
# This is a dummy build to get the dependencies cached.
RUN cargo build --release
#
# Now copy in the rest of the sources.
COPY build.rs Cargo.toml Cargo.lock /usr/src/varfish-server-worker/
COPY src /usr/src/varfish-server-worker/src/
COPY protos /usr/src/varfish-server-worker/protos/
# COPY utils/alpine-linker-script.sh /usr/src/varfish-server-worker/utils/
# RUN chmod a+rx /usr/src/varfish-server-worker/utils/alpine-linker-script.sh
# COPY .cargo /usr/src/varfish-server-worker/.cargo/
## Touch main.rs to prevent cached release build.
RUN touch /usr/src/varfish-server-worker/src/main.rs
# This is the actual application build.
RUN cargo build --release

# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------

FROM debian:bookworm-slim AS runtime

# Install dependencies (and cleanup afterwards)
RUN apt-get update && \
    apt-get install -y libsqlite3-0 && \
    apt-get clean autoclean && \
    apt-get autoremove --yes && \
    rm -rf /var/lib/{apt,dpkg,cache,log}

# Copy application binary from builder image
COPY --from=builder \
    /usr/src/varfish-server-worker/target/release/varfish-server-worker \
    /usr/local/bin

# Copy the entrypoint script and make it executable.
COPY utils/docker/entrypoint.sh /
RUN chmod a+rx /entrypoint.sh

# Set the entrypoint.
ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]
