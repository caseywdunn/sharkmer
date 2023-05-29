# If data/SRR5324768_pass_1.fastq does not exist, uncompress data/SRR5324768_pass_1.fastq.gz
if [ ! -f data/SRR5324768_pass_1.fastq ]; then
	gunzip -c data/SRR5324768_pass_1.fastq.gz > data/SRR5324768_pass_1.fastq
fi

# Build 
cargo build --release
mkdir -p $PREFIX/bin
cp target/release/sharkmer $PREFIX/bin