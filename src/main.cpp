#include <iostream>
#include <cstdint>
#include <cassert>
#include <sstream>
#include <vector>
#include <filesystem>

#include "3rdparty/CLI11.hpp"
#include "3rdparty/spdlog/spdlog.h"
#include "3rdparty/BS_thread_pool.hpp"

#include "mceliece.h"
#include "candidates.h"
#include "leakfile.h"
#include "poly.h"
#include "classic.h"
#include "botan.h"

#include "main.h"

struct CommandSetup {
    std::string kem;

    CommandSetup() {
    }

    template<McElieceParameters P>
    void op_setup_classic() {
        McEliece<P> params;
        Candidates<decltype(params)> classic_candidates;
    }
    template<McElieceParameters P>
    void op_setup_botan() {
        McEliece<P> params;
        Candidates<decltype(params)> botan_candidates(true);
    }

    void run() {
        if (kem == "mceliece348864") {
            op_setup_classic<Classic348864>();
        } else if (kem == "mceliece460896") {
            op_setup_classic<Classic460896>();
        } else if (kem == "mceliece6688128") {
            op_setup_classic<Classic6688128>();
        } else if (kem == "mceliece6960119") {
            op_setup_classic<Classic6960119>();
        } else if (kem == "mceliece8192128") {
            op_setup_classic<Classic8192128>();
        } else if (kem == "mceliece51220") {
            op_setup_classic<Classic51220>();
        } else if (kem == "mceliece102450") {
            op_setup_classic<Classic102450>();
        } else if (kem == "botan51220") {
            op_setup_botan<Botan51220>();
        } else if (kem == "botan102450") {
            op_setup_botan<Botan102450>();
        } else if (kem == "botan348864") {
            op_setup_botan<Botan348864>();
        } else if (kem == "botan460896") {
            op_setup_botan<Botan460896>();
        } else if (kem == "botan6688128") {
            op_setup_botan<Botan6688128>();
        } else if (kem == "botan6960119") {
            op_setup_botan<Botan6960119>();
        } else if (kem == "botan8192128") {
            op_setup_botan<Botan8192128>();
        }
    }
};

struct CommandRecover {
    std::string filename;
    bool verbose;
    unsigned int seed;
    double readout_error;
    bool with_distances;
    std::string distances_fname;
    bool with_report;
    std::string report_fname;
    bool external_error;

    CommandRecover() {
        verbose = false;
        seed = 0;
        readout_error = 0.3;
        with_distances = false;
        distances_fname = "distances-report.csv";
        with_report = false;
        report_fname = "report.json";
	external_error = false;
    }

    template<McElieceParameters P>
    void op_recover_classic(Leakfile& leakfile) {
        McEliece<P> params;
        Classic<decltype(params)> classic;
        classic.verbose = verbose;
        classic.recover_from_leak(leakfile, seed, readout_error, external_error, with_distances, distances_fname,
                                  with_report, report_fname);
    }
    template<McElieceParameters P>
    void op_recover_botan(Leakfile& leakfile) {
        McEliece<P> params;
        Botan<decltype(params)> botan;
        botan.verbose = verbose;
	if (external_error) {
		spdlog::critical("external error not supported for botan key recovery.");
		return;
	}
        botan.recover_from_leak(leakfile, seed, readout_error, with_report, report_fname);
    }

    void run() {
        Leakfile leakfile(filename);
        if (leakfile.kem == "mceliece348864") {
            op_recover_classic<Classic348864>(leakfile);
        } else if (leakfile.kem == "mceliece460896") {
            op_recover_classic<Classic460896>(leakfile);
        } else if (leakfile.kem == "mceliece6688128") {
            op_recover_classic<Classic6688128>(leakfile);
        } else if (leakfile.kem == "mceliece6960119") {
            op_recover_classic<Classic6960119>(leakfile);
        } else if (leakfile.kem == "mceliece8192128") {
            op_recover_classic<Classic8192128>(leakfile);
        } else if (leakfile.kem == "mceliece51220") {
            op_recover_classic<Classic51220>(leakfile);
        } else if (leakfile.kem == "mceliece102450") {
            op_recover_classic<Classic102450>(leakfile);
        } else if (leakfile.kem == "botan51220") {
            op_recover_botan<Botan51220>(leakfile);
        } else if (leakfile.kem == "botan102450") {
            op_recover_botan<Botan102450>(leakfile);
        } else if (leakfile.kem == "botan348864") {
            op_recover_botan<Botan348864>(leakfile);
        } else if (leakfile.kem == "botan460896") {
            op_recover_botan<Botan460896>(leakfile);
        } else if (leakfile.kem == "botan6688128") {
            op_recover_botan<Botan6688128>(leakfile);
        } else if (leakfile.kem == "botan6960119") {
            op_recover_botan<Botan6960119>(leakfile);
        } else if (leakfile.kem == "botan8192128") {
            op_recover_botan<Botan8192128>(leakfile);
        }
    }
};

BS::thread_pool pool;

int main(int argc, char* argv[]) {
    spdlog::info("Welcome to Stockfish!");
    CLI::App app{"Artefacts for the McEliece paper"};
    bool debug = false;

    app.add_flag_callback("--debug", []() { spdlog::set_level(spdlog::level::debug); }, "more verbose output")->capture_default_str();
    app.add_option_function<int>("--threads", [](const int&nr_threads) { pool.reset(nr_threads); }, "number of parallel threads")->default_val(pool.get_thread_count());

    auto cmd_recover = std::make_shared<CommandRecover>();
    auto* sub_recover = app.add_subcommand("recover", "calculate secret from leakfile");
    sub_recover->add_option("--seed", cmd_recover->seed, "random seed (0 means OS random)")->capture_default_str();
    sub_recover->add_option("--error", cmd_recover->readout_error, "readout error rate")->capture_default_str();
    sub_recover->add_flag("--external-error", cmd_recover->external_error, "error is externally applied")->capture_default_str();
    sub_recover->add_flag("--with-distances", cmd_recover->with_distances, "generate distances report")->capture_default_str();
    sub_recover->add_option("--distances-file", cmd_recover->distances_fname, "generate distances report in file")->capture_default_str();
    sub_recover->add_flag("--with-report", cmd_recover->with_report, "generate result report in file")->capture_default_str();
    sub_recover->add_option("--report-file", cmd_recover->report_fname, "if result is generated, store it in file")->capture_default_str();
    sub_recover->add_option("LEAKFILE", cmd_recover->filename, "filename of leakfile")->required();
    // Set the run function as callback to be called when this subcommand is issued.
    sub_recover->callback([cmd_recover]() { cmd_recover->run(); });

    auto cmd_setup = std::make_shared<CommandSetup>();
    auto* sub_setup = app.add_subcommand("setup", "pregenerate cache file");
    sub_setup->add_option("KEM", cmd_setup->kem, "kem name")->required();
    // Set the run function as callback to be called when this subcommand is issued.
    sub_setup->callback([cmd_setup]() { cmd_setup->run(); });

    app.parse_complete_callback([](){
        spdlog::info("Using {0} threads for parallel processing.", pool.get_thread_count());
    });

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    return 0;
}
