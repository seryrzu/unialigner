build_type?="Release"

.PHONY: clean cmake all

all: tandem_aligner

cmake:
	mkdir -p build
	cd build && cmake .. -DCMAKE_BUILD_TYPE="${build_type}"

tandem_aligner: cmake
	$(MAKE) -C build all
	mkdir -p build/bin
	mv $(abspath build/src/projects/tandem_aligner/tandem_aligner) build/bin/tandem_aligner

test_launch: tandem_aligner
	build/bin/tandem_aligner \
		--first test_dataset/first.fasta \
		--second test_dataset/second.fasta \
		-o test_dataset/test_launch \
		--debug
	grep -q "Thank you for using TandemAligner!" test_dataset/test_launch/tandem_aligner.log
	diff test_dataset/test_launch/cigar.txt test_dataset/true_cigar.txt

clean:
	-rm -r build
	-rm -r test_dataset/test_launch