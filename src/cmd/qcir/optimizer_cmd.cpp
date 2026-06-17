/****************************************************************************
  PackageName  [ qcir/optimizer ]
  Synopsis     [ Define optimizer package commands ]
  Author       [ Design Verification Lab ]
  Copyright    [ Copyright(c) 2023 DVLab, GIEE, NTU, Taiwan ]
****************************************************************************/

#include <spdlog/spdlog.h>

#include <sstream>
#include <string>

#include "argparse/arg_type.hpp"
#include "cli/cli.hpp"
#include "cmd/qcir_mgr.hpp"
#include "qcir/optimizer/optimizer.hpp"
#include "qcir/qcir.hpp"
#include "tableau/phasepoly/config.hpp"
#include "tableau/phasepoly/qcir_optimizer.hpp"
#include "util/data_structure_manager_common_cmd.hpp"
#include "util/dvlab_string.hpp"
#include "util/phase.hpp"

using namespace dvlab::argparse;
using dvlab::CmdExecResult;
using dvlab::Command;

extern bool stop_requested();

namespace qsyn::qcir {

//----------------------------------------------------------------------
//    Optimize
//----------------------------------------------------------------------
Command qcir_optimize_cmd(QCirMgr& qcir_mgr) {
    return {"optimize",
            [](ArgumentParser& parser) {
                parser.description("optimize QCir");

                parser.add_argument<std::string>("strategy")
                    .help("optimization strategy")
                    .default_value("basic")
                    .constraint(choices_allow_prefix(
                        {"basic",
                         "teleport",
                         "blaqsmith",
                         "phasepoly"}));

                parser.add_argument<double>("--init-temp")
                    .default_value(0.5)
                    .help("initial temperature for annealing");

                parser.add_argument<bool>("-p", "--physical")
                    .default_value(false)
                    .action(store_true)
                    .help("optimize physical circuit, i.e preserve the swap path");
                parser.add_argument<bool>("-c", "--copy")
                    .default_value(false)
                    .action(store_true)
                    .help("copy a circuit to perform optimization");
                parser.add_argument<bool>("-s", "--statistics")
                    .default_value(false)
                    .action(store_true)
                    .help("count the number of rules operated in optimizer.");
                parser.add_argument<bool>("-t", "--tech")
                    .default_value(false)
                    .action(store_true)
                    .help("Only perform optimizations preserving gate sets and qubit connectivities.");

                parser.add_argument<size_t>("--max-queue")
                    .default_value(10000)
                    .help("PhasePoly: cap on the A* open set (default: 10000).");
                parser.add_argument<size_t>("--max-expansions")
                    .default_value(500000)
                    .help("PhasePoly: hard cap on state expansions (default: 500000).");
                parser.add_argument<std::string>("--group-size")
                    .default_value("1,2,3,5")
                    .help("PhasePoly: comma-separated adjacent block group sizes to try (default: 1,2,3,5).");
            },
            [&](ArgumentParser const& parser) {
                if (!dvlab::utils::mgr_has_data(qcir_mgr)) return CmdExecResult::error;
                Optimizer optimizer;
                std::optional<QCir> result;
                std::string procedure_str{};

                enum class Strategy {
                    basic,
                    teleport,
                    blaqsmith,
                    phasepoly
                };

                auto strategy = [&]() -> Strategy {
                    if (dvlab::str::is_prefix_of(parser.get<std::string>("strategy"), "teleport")) {
                        return Strategy::teleport;
                    }
                    if (dvlab::str::is_prefix_of(parser.get<std::string>("strategy"), "blaqsmith")) {
                        return Strategy::blaqsmith;
                    }
                    if (dvlab::str::is_prefix_of(parser.get<std::string>("strategy"), "phasepoly")) {
                        return Strategy::phasepoly;
                    }
                    return Strategy::basic;
                }();

                switch (strategy) {
                    case Strategy::teleport:
                        phase_teleport(*qcir_mgr.get());
                        procedure_str = "Phase Teleport";
                        break;
                    case Strategy::blaqsmith:
                        optimize_2q_count(*qcir_mgr.get(), parser.get<double>("--init-temp"), 2, 2);
                        procedure_str = "Blaqsmith";
                        break;
                    case Strategy::phasepoly: {
                        experimental::phasepoly::PhasePolyConfig config;
                        config.max_queue_size  = parser.get<size_t>("--max-queue");
                        config.max_expansions  = parser.get<size_t>("--max-expansions");

                        config.group_sizes.clear();
                        {
                            std::stringstream ss(parser.get<std::string>("--group-size"));
                            std::string token;
                            while (std::getline(ss, token, ',')) {
                                if (!token.empty()) {
                                    config.group_sizes.push_back(std::stoull(token));
                                }
                            }
                        }
                        if (config.group_sizes.empty()) {
                            config.group_sizes = {1, 2, 3, 5};
                        }

                        experimental::phasepoly::PhasePolyOptimizeStats stats{};
                        result = experimental::phasepoly::optimize_qcir_phasepoly(
                            *qcir_mgr.get(), config, true, &stats);

                        if (!result) {
                            spdlog::error("PhasePoly optimization failed.");
                            return CmdExecResult::error;
                        }

                        spdlog::info(
                            "PhasePoly: {} block(s), {} improved, CX {} -> {}",
                            stats.num_blocks,
                            stats.blocks_optimized,
                            stats.cx_before,
                            stats.cx_after);
                        procedure_str = "PhasePolyOpt";
                        if (parser.get<bool>("--copy")) {
                            qcir_mgr.add(qcir_mgr.get_next_id(), std::make_unique<QCir>(std::move(*result)));
                        } else {
                            qcir_mgr.set(std::make_unique<QCir>(std::move(*result)));
                        }
                        break;
                    }
                    case Strategy::basic: {
                        if (parser.get<bool>("--tech") || !qcir_mgr.get()->get_gate_set().empty()) {
                            result        = optimizer.trivial_optimization(*qcir_mgr.get());
                            procedure_str = "Tech Optimize";
                        } else {
                            result        = optimizer.basic_optimization(*qcir_mgr.get(), {.doSwap          = !parser.get<bool>("--physical"),
                                                                                           .maxIter         = 1000,
                                                                                           .printStatistics = parser.get<bool>("--statistics")});
                            procedure_str = "Optimize";
                        }
                        if (result == std::nullopt) {
                            spdlog::error("Fail to optimize circuit.");
                            return CmdExecResult::error;
                        }

                        if (parser.get<bool>("--copy")) {
                            qcir_mgr.add(qcir_mgr.get_next_id(), std::make_unique<QCir>(std::move(*result)));
                        } else {
                            qcir_mgr.set(std::make_unique<QCir>(std::move(*result)));
                        }
                    }
                }
                if (stop_requested()) {
                    procedure_str += "[INT]";
                }
                qcir_mgr.get()->add_procedure(procedure_str);

                return CmdExecResult::done;
            }};
}

}  // namespace qsyn::qcir
