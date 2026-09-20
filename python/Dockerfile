ARG RUST_VERSION=1.98.0
ARG RUST_TOOLCHAIN=nightly-2026-06-26
FROM rust:${RUST_VERSION}-bookworm

ARG RUST_TOOLCHAIN
RUN rustup toolchain install "${RUST_TOOLCHAIN}" --profile minimal \
    --component rustfmt --component clippy
ENV RUSTUP_TOOLCHAIN=${RUST_TOOLCHAIN}

ENV CARGO_TERM_COLOR=always \
    PYO3_PYTHON=/usr/bin/python3 \
    VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH

RUN apt-get update \
    && apt-get install --no-install-recommends -y python3 python3-venv python3-dev pkg-config libssl-dev zlib1g-dev \
    && python3 -m venv /opt/venv \
    && pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir maturin==1.11.5 pytest==8.4.2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
COPY Cargo.toml Cargo.lock ./
COPY python/Cargo.toml python/Cargo.lock python/pyproject.toml python/PYPI-README.md python/
RUN mkdir src python/src benches \
    && printf 'pub fn dependency_cache_probe() {}\n' > src/lib.rs \
    && for bench in build-performance compression query-performance; do \
        printf 'fn main() {}\n' > "benches/${bench}.rs"; \
    done \
    && printf 'pub fn dependency_cache_probe() {}\n' > python/src/lib.rs \
    && cargo build --manifest-path python/Cargo.toml --release --locked --features abi3,extension-module \
    && rm -rf src python/src benches

COPY . .
RUN touch src/lib.rs src/main.rs python/src/lib.rs \
    && maturin build --manifest-path python/Cargo.toml --release --locked \
        --features abi3,extension-module --compatibility pypi --out /tmp/wheels \
    && pip install --no-cache-dir /tmp/wheels/*.whl \
    && pytest -q python/tests \
    && useradd --create-home --shell /bin/bash gai \
    && chown -R gai:gai /workspace /opt/venv

USER gai
CMD ["pytest", "-q", "python/tests"]
