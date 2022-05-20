//
// Created by Andrey Bzikadze on 05/05/22.
//

#pragma once

#include <common/logging.hpp>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include "lcp_interval.hpp"
#include "min_interval.hpp"
#include "sparse_aligner.hpp"

namespace minseq {

class MinSeqAligner {
    logging::Logger &logger;
    int max_freq{1};
    const std::experimental::filesystem::path output_dir;

    [[nodiscard]] std::string ReadContig(const std::experimental::filesystem::path &path) const {
        io::SeqReader reader(path);
        std::vector<Contig> vec{reader.readAllContigs()};
        VERIFY(vec.size()==1);
        Contig contig{std::move(vec[0])};
        logger.info() << "Length " << contig.seq.size() << ", name "
                      << contig.id
                      << " from " << path << std::endl;
        return contig.str();
    }

    [[nodiscard]] std::string ConcatContigs(const std::string &first,
                                            const std::string &second) const {
        std::stringstream concat_stream;
        concat_stream << first << '$' << second << '#';
        return concat_stream.str();
    }

 public:
    MinSeqAligner(logging::Logger &logger,
                  std::experimental::filesystem::path output_dir,
                  const int max_freq) :
        logger{logger}, output_dir{std::move(output_dir)}, max_freq{max_freq} {}

    void Find(const std::experimental::filesystem::path &first_path,
              const std::experimental::filesystem::path &second_path) const {
        io::SeqReader first_reader(first_path);
        std::vector<Contig> first_vec{first_reader.readAllContigs()};
        VERIFY(first_vec.size()==1);
        const std::string first = ReadContig(first_path);
        const std::string second = ReadContig(second_path);
        const std::string concat = ConcatContigs(first, second);

        logger.info() << "Building suffix array...\n";
        const suffix_array::SuffixArray<std::string> suf_arr(concat);
        logger.info() << "Building LCP array...\n";
        const suffix_array::LCP<std::string> lcp(suf_arr);

        MinIntervalFinder segment_finder(max_freq);
        logger.info() << "Computing rare segments...\n";
        const MinIntervalCollections
            int_col = segment_finder.Find(lcp, first.size());
        std::ofstream os(output_dir/"shortest_matches.tsv");
        os << int_col;

        logger.info() << "Aligning...\n";
        const Cigar cigar = SparseAligner(logger).Align(int_col, first, second);
        std::ofstream cigar_os(output_dir/"cigar.txt");
        cigar_os << cigar;
        logger.info() << "Cigar exported\n";
    }
};

} // namespace minseq

