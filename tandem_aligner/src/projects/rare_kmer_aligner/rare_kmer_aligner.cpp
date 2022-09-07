#include <iomanip>
#include <iostream>

#include <common/cl_parser.hpp>
#include <common/logging.hpp>

#include "rare_kmer_aligner.hpp"

int main(int argc, char ** argv) {
    CLParser parser {{"output-dir=", "first=", "second=",
                      "k-mer-size=10", "max_rare_freq=1", "tol_gap=100"}, {},
                     {"o=output-dir", "f=first", "s=second",
                      "k=k-mer-size", "r=max_rare_freq", "g=tol_gap"}};
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cerr << "Incorrect parameters" << std::endl;
        std::cerr << parser.check() << std::endl;
        return 1;
    }

    const std::experimental::filesystem::path output_dir{parser.getValue("output-dir")};
    ensure_dir_existance(output_dir);
    logging::LoggerStorage ls{output_dir, "rare_kmer_aligner"};
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());

    auto time_point {std::chrono::system_clock::now()};
    std::time_t now = std::chrono::system_clock::to_time_t(time_point);
    logger << "Time: " << std::put_time(std::localtime(&now), "%c %Z") << std::endl;

    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;

    const std::experimental::filesystem::path first_path =
            std::experimental::filesystem::canonical(parser.getValue("first"));
    const std::experimental::filesystem::path second_path =
            std::experimental::filesystem::canonical(parser.getValue("second"));

    size_t k = std::stoi(parser.getValue("k-mer-size"));
    size_t max_rare_freq = std::stoi(parser.getValue("max_rare_freq"));
    size_t tol_gap = std::stoi(parser.getValue("tol_gap"));

    rare_kmer_aligner::rare_kmer_align<unsigned __int128>(first_path,
                                                          second_path,
                                                          k,
                                                          max_rare_freq,
                                                          tol_gap,
                                                          output_dir,
                                                          logger);

    return 0;
}
