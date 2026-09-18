FROM rust:1.95-bookworm

WORKDIR /app

COPY Cargo.lock Cargo.toml ./
COPY src ./src

RUN cargo build --release
