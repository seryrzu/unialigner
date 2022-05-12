//
// Created by Andrey Bzikadze on 05/05/22.
//

#include <iomanip>
#include <iostream>

#include <common/cl_parser.hpp>
#include <common/logging.hpp>

#include "dotplot.hpp"

int main(int argc, char **argv) {
    CLParser parser{{"output-dir=", "first=", "second="}, {},
                    {"o=output-dir", "f=first", "s=second"}};
    parser.parseCL(argc, argv);
    if (!parser.check().empty()) {
        std::cerr << "Incorrect parameters" << std::endl;
        std::cerr << parser.check() << std::endl;
        return 1;
    }

    const std::experimental::filesystem::path
        output_dir{parser.getValue("output-dir")};
    ensure_dir_existance(output_dir);
    logging::LoggerStorage ls{output_dir, "tandem_aligner"};
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile());

    auto time_point{std::chrono::system_clock::now()};
    std::time_t now = std::chrono::system_clock::to_time_t(time_point);
    logger << "Time: " << std::put_time(std::localtime(&now), "%c %Z")
           << std::endl;

    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;

    const std::experimental::filesystem::path first_path =
        std::experimental::filesystem::canonical(parser.getValue("first"));
    const std::experimental::filesystem::path second_path =
        std::experimental::filesystem::canonical(parser.getValue("second"));

    dotplot::DotPlotConstructor(logger, output_dir)
        .Construct(first_path, second_path);

    logger.info() << "Thank you for using TandemAligner!\n";
    return 0;
}
