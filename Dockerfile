ARG RUST_VERSION=1.98.0
ARG RUST_TOOLCHAIN=nightly-2026-06-26
FROM rust:${RUST_VERSION}-bookworm

ARG RUST_TOOLCHAIN
RUN rustup toolchain install "${RUST_TOOLCHAIN}" --profile minimal \
    --component rustfmt --component clippy
ENV RUSTUP_TOOLCHAIN=${RUST_TOOLCHAIN}

ENV CARGO_TERM_COLOR=always \
    RUST_BACKTRACE=1

RUN apt-get update \
    && apt-get install --no-install-recommends -y pkg-config libssl-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY Cargo.lock Cargo.toml ./
RUN mkdir src benches \
    && printf 'pub fn dependency_cache_probe() {}\n' > src/lib.rs \
    && for bench in build-performance compression query-performance; do \
        printf 'fn main() {}\n' > "benches/${bench}.rs"; \
    done \
    && cargo build --release --locked --lib \
    && rm -rf src benches

COPY . .

RUN touch src/lib.rs src/main.rs \
    && cargo build --release --locked \
    && cargo test --locked --all-targets --all-features

RUN useradd --create-home --shell /bin/bash gni \
    && chown -R gni:gni /app
USER gni

CMD ["bash"]
