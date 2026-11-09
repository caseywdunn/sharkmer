# Source this file to use the project-local Rust and Cargo installation.
export CARGO_HOME=/nfs/roberts/project/pi_cwd7/cwd7/sharkmer/tmp/uncertain-products-20260927/toolchain/cargo
export RUSTUP_HOME=/nfs/roberts/project/pi_cwd7/cwd7/sharkmer/tmp/uncertain-products-20260927/toolchain/rustup
case ":${PATH}:" in
    *":${CARGO_HOME}/bin:"*) ;;
    *) export PATH="${CARGO_HOME}/bin:${PATH}" ;;
esac
