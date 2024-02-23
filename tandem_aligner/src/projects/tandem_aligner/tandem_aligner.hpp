//
// Created by Andrey Bzikadze on 05/05/22.
//

#pragma once

#include <queue>
#include <common/logging.hpp>
#include <sequences/contigs.hpp>
#include <sequences/seqio.hpp>
#include "lcp_interval.hpp"
#include "min_interval.hpp"
#include "sparse_aligner.hpp"

namespace tandem_aligner {

struct MinSeqTask {
    std::list<CigarFragment>::iterator cigar_it;
    int64_t st1{0}, len1{0};
    int64_t st2{0}, len2{0};
};

class TandemAligner {
    logging::Logger &logger;
    int max_freq{1};
    bool force_highfreq_search{false};
    const std::experimental::filesystem::path output_dir;
    bool no_paths;
    bool bridges;

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

    void RunTask(std::queue<MinSeqTask> &queue,
                 Cigar &main_cigar,
                 const std::string &first,
                 const std::string &second,
                 bool exprt,
                 bool assert_validity = true) const {
        MinSeqTask task = queue.front();
        const std::string first_substr = first.substr(task.st1, task.len1);
        const std::string
            second_substr = second.substr(task.st2, task.len2);
        const std::string
            concat = ConcatContigs(first_substr, second_substr);

        logger.debug() << "Building suffix array...\n";
        const suffix_array::SuffixArray<std::string> suf_arr(concat);
        logger.debug() << "Building LCP array...\n";
        const suffix_array::LCP<std::string> lcp(suf_arr);

        MaxDisjointIntervalFinder
            segment_finder(max_freq, force_highfreq_search, exprt,
                           output_dir/"min_interval_finder");
        logger.debug() << "Computing rare segments...\n";
        const MaxDisjointIntervalCollections
            int_col = segment_finder.Find(lcp, task.len1);

        if (exprt) {
            std::ofstream os(output_dir/"shortest_matches.tsv");
            os << int_col;
        }

        logger.debug() << "Aligning...\n";
        Cigar cigar = SparseAligner(logger, output_dir, exprt & no_paths, exprt & bridges).Align(int_col,
                                                  first_substr, second_substr); //exprt & for primary alignment
        auto main_cigar_it = main_cigar.AssignInterval(cigar, task.cigar_it);
        logger.debug() << "Finished alignment\n";
        if (cigar.Size() > 2) {
            int64_t i{task.st1}, j{task.st2};
            auto it1 = cigar.cbegin(), it2 = ++cigar.cbegin();
            for (; it2!=cigar.cend(); ++it1, ++it2, ++main_cigar_it) {
                const CigarFragment &fragment1 = *it1;
                const CigarFragment &fragment2 = *it2;
                if ((fragment1.mode==CigarMode::D
                    or fragment1.mode==CigarMode::I) and
                    (fragment2.mode==CigarMode::D
                        or fragment2.mode==CigarMode::I) and
                    std::min(fragment1.length, fragment2.length) > 1) {
                    if (fragment1.mode==CigarMode::D) {
                        queue.push({main_cigar_it,
                                    i, int64_t(fragment1.length),
                                    j, int64_t(fragment2.length)});
                    } else {
                        queue.push({main_cigar_it,
                                    i, int64_t(fragment2.length),
                                    j, int64_t(fragment1.length)});
                    }
                }
                if (fragment1.mode==CigarMode::M
                    or fragment1.mode==CigarMode::X) {
                    i += fragment1.length;
                    j += fragment1.length;
                    continue;
                } else if (fragment1.mode==CigarMode::D) {
                    i += fragment1.length;
                } else {
                    j += fragment1.length;
                }
            }
        }
        if (assert_validity)
            cigar.AssertValidity(first, second);
    }

    void AssignMismatches(Cigar &cigar,
                          const std::string &first,
                          const std::string &second) const {
        if (cigar.Size() < 2)
            return;
        auto it1 = cigar.begin(), it2 = ++cigar.begin();
        int i{0}, j{0};
        std::unordered_map<int, int> counter;
        for (; it2!=cigar.end(); ++it1, ++it2) {
            int64_t length1{it1->length}, length2{it2->length};
            CigarMode mode1{it1->mode}, mode2{it2->mode};
            if ((mode1==CigarMode::D or mode1==CigarMode::I) and
                (mode2==CigarMode::D or mode2==CigarMode::I) and
                length1==length2) {
                counter[length1]++;
                it2 = cigar.Erase(it2);
                it1 = cigar.Erase(it1);
                VERIFY(it1==it2);

                int64_t run_len{1}, k{1};
                bool is_eq{first[i]==second[j]};
                i++, j++;
                while (true) {
                    if (is_eq!=(first[i]==second[j]) or k==length1) {
                        it2 = cigar.Insert(it2,
                                           {run_len, is_eq ? CigarMode::M
                                                           : CigarMode::X});
                        it1 = it2++;
                        if (k==length1)
                            break;
                        is_eq = !is_eq;
                        run_len = 0;
                    }
                    run_len++, k++, i++, j++;
                }
            } else if (mode1==CigarMode::M or mode1==CigarMode::X) {
                i += length1;
                j += length1;
            } else if (mode1==CigarMode::D) {
                i += length1;
            } else {
                j += length1;
            }
        }
        for (auto [length, cnt] : counter) {
            logger.trace() << "Square Indel-block of length " << length
                << " appears " << cnt << " times\n";
        }
        cigar.AssertValidity(first, second);
    }

 public:
    TandemAligner(logging::Logger &logger,
                  std::experimental::filesystem::path output_dir,
                  const int max_freq,
                  const bool force_highfreq_search,
                  bool no_paths,
                  bool bridges) :
        logger{logger}, output_dir{std::move(output_dir)}, max_freq{max_freq},
        force_highfreq_search{force_highfreq_search}, no_paths{no_paths}, bridges{bridges} {}

    void Find(const std::experimental::filesystem::path &first_path,
              const std::experimental::filesystem::path &second_path) const {
        io::SeqReader first_reader(first_path);
        std::vector<Contig> first_vec{first_reader.readAllContigs()};
        VERIFY(first_vec.size()==1);
        const std::string first = ReadContig(first_path);
        const std::string second = ReadContig(second_path);

        Cigar cigar;
        std::queue<MinSeqTask> queue;
        bool matches_exported{false};
        queue.push({cigar.begin(), 0, (int64_t) first.size(), 0,
                    (int64_t) second.size()});
        logger.info() << "Running primary alignment...\n";
        if (no_paths) std::ofstream no_paths(output_dir/"no_paths.csv"); //just empty the file 
        if (bridges) std::ofstream no_paths(output_dir/"bridges.txt"); //just empty the file 
        RunTask(queue, cigar, first, second,
                /*export_matches*/ true);
        logger.info() << "Number of indel-blocks " << queue.size() << "\n";
        queue.pop();
        logger.info() << "Finished running primary alignment\n";
        cigar.Summary(logger);
        logger.info() << "Running recursive alignments...\n";

        {
            std::string cigar_outfile = output_dir/"cigar_primary.txt";
            std::ofstream cigar_os(cigar_outfile);
            cigar_os << cigar;
            logger.info() << "Primary cigar exported to " << cigar_outfile
                          << "\n\n";
        }

        for (; not queue.empty(); queue.pop()) {
            RunTask(queue, cigar, first, second,
                /*export_matches*/ false,
                /*assert_validity*/ false);
        }
        cigar.AssertValidity(first, second);
        logger.info() << "Finished running recursive alignment\n";
        cigar.Summary(logger);

        {
            std::string cigar_outfile = output_dir/"cigar_recursive.txt";
            std::ofstream cigar_os(cigar_outfile);
            cigar_os << cigar;
            logger.info() << "Cigar after recursion exported to "
                          << cigar_outfile << "\n\n";
        }

        logger.info() << "Assigning mismatches...\n";
        AssignMismatches(cigar, first, second);
        logger.info() << "Finished assigning mismatches\n";
        cigar.Summary(logger);

        {
            std::string cigar_outfile = output_dir/"cigar.txt";
            std::ofstream cigar_os(cigar_outfile);
            cigar_os << cigar;
            logger.info() << "Cigar w/ mismatches exported to " << cigar_outfile
                          << "\n";
        }
    }
};

} // namespace tandem_aligner
