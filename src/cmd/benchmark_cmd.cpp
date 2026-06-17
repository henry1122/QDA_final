/****************************************************************************
  PackageName  [ cmd ]
  Synopsis     [ PhasePoly Table 1 reproduction benchmark ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2025 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include "./benchmark_cmd.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "argparse/arg_type.hpp"
#include "cli/cli.hpp"
#include "tableau/phasepoly/config.hpp"
#include "tableau/phasepoly/table1_benchmark.hpp"

using namespace dvlab::argparse;
using dvlab::CmdExecResult;
using dvlab::Command;

namespace qsyn::experimental::phasepoly {

namespace {

Command benchmark_table1_cmd() {
    return {"table1",
            [](ArgumentParser& parser) {
                parser.description(
                    "Reproduce PhasePoly paper Table 1: Todd (+P) then "
                    "synthesize with pp / mst / gray / gstair / naive.");

                parser.add_argument<std::string>("--bench-dir")
                    .default_value("benchmark/qc/optimized")
                    .help("directory containing *_pyzx.qc circuits");

                parser.add_argument<std::string>("--format")
                    .default_value("markdown")
                    .constraint(choices_allow_prefix({"markdown", "csv", "phase-csv"}))
                    .help("output format: markdown, csv (full circuit), phase-csv (CX+Rz blocks only)");

                parser.add_argument<std::string>("--output")
                    .default_value("")
                    .help("write results to file (default: stdout)");

                parser.add_argument<std::string>("--phase-output")
                    .default_value("")
                    .help("also write phase-only CX+Rz CSV to this path");

                parser.add_argument<std::string>("--circuits")
                    .default_value("")
                    .help("comma-separated circuit filenames (default: full Table 1 set)");

                parser.add_argument<size_t>("--max-queue")
                    .default_value(10000)
                    .help("PhasePoly A* queue cap (default: 10000)");

                parser.add_argument<size_t>("--max-expansions")
                    .default_value(500000)
                    .help("PhasePoly expansion cap (default: 500000)");

                parser.add_argument<std::string>("--group-size")
                    .default_value("1,2,3,5")
                    .help("multiblock group sizes for pp column (default: 1,2,3,5)");

                parser.add_argument<bool>("--no-todd")
                    .default_value(false)
                    .help("disable block-level Todd preprocessing (enabled by default)");
            },
            [&](ArgumentParser const& parser) {
                PhasePolyConfig config;
                config.max_queue_size = parser.get<size_t>("--max-queue");
                config.max_expansions = parser.get<size_t>("--max-expansions");

                Table1BenchmarkOptions options;
                options.use_todd = !parser.get<bool>("--no-todd");

                config.group_sizes.clear();
                {
                    std::stringstream ss(parser.get<std::string>("--group-size"));
                    std::string token;
                    while (std::getline(ss, token, ',')) {
                        if (!token.empty()) config.group_sizes.push_back(std::stoull(token));
                    }
                }
                if (config.group_sizes.empty()) config.group_sizes = {1, 2, 3, 5};

                std::vector<std::filesystem::path> paths;
                auto const bench_dir = std::filesystem::path{parser.get<std::string>("--bench-dir")};

                auto const names = [&]() -> std::vector<std::string> {
                    auto const list = parser.get<std::string>("--circuits");
                    if (list.empty()) return table1_circuit_names();
                    std::vector<std::string> out;
                    std::stringstream ss(list);
                    std::string token;
                    while (std::getline(ss, token, ',')) {
                        if (!token.empty()) out.push_back(token);
                    }
                    return out;
                }();

                for (auto const& name : names) {
                    paths.push_back(bench_dir / name);
                }

                auto const summary = benchmark_table1(paths, config, options);

                std::ostringstream buffer;
                if (parser.get<std::string>("--format") == "csv") {
                    print_table1_csv(summary, buffer);
                } else if (parser.get<std::string>("--format") == "phase-csv") {
                    print_table1_phase_csv(summary, buffer);
                } else {
                    print_table1_markdown(summary, buffer);
                }

                auto const out_path = parser.get<std::string>("--output");
                if (out_path.empty()) {
                    fmt::print("{}", buffer.str());
                } else {
                    std::ofstream ofs{out_path};
                    if (!ofs) {
                        spdlog::error("Cannot write {}", out_path);
                        return CmdExecResult::error;
                    }
                    ofs << buffer.str();
                    spdlog::info("Table 1 results written to {}", out_path);
                }

                auto const phase_out_path = parser.get<std::string>("--phase-output");
                if (!phase_out_path.empty()) {
                    std::ostringstream phase_buffer;
                    print_table1_phase_csv(summary, phase_buffer);
                    std::ofstream pofs{phase_out_path};
                    if (!pofs) {
                        spdlog::error("Cannot write {}", phase_out_path);
                        return CmdExecResult::error;
                    }
                    pofs << phase_buffer.str();
                    spdlog::info("Table 1 phase-only results written to {}", phase_out_path);
                }

                return CmdExecResult::done;
            }};
}

}  // namespace

bool add_benchmark_cmds(dvlab::CommandLineInterface& cli) {
    auto cmd = Command{
        "benchmark",
        [](ArgumentParser& parser) {
            parser.description("PhasePoly paper benchmarks");
            parser.add_subparsers("subcmd").required(true);
        },
        [](ArgumentParser const& /*parser*/) { return CmdExecResult::error; }};

    cmd.add_subcommand("subcmd", benchmark_table1_cmd());
    return cli.add_command(cmd);
}

}  // namespace qsyn::experimental::phasepoly
