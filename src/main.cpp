#include <volt/cli/common.h>
#include <volt/coordination_service.h>

using namespace Volt;
using namespace Volt::CLI;

void showUsage(const std::string& name) {
    printUsageHeader(name, "Volt - Coordination Analysis");
    std::cerr
        << "  --cutoff <float>              Cutoff radius for neighbor search. [default: 3.2]\n"
        << "  --rdfBins <int>               Number of bins for RDF calculation. [default: 500]\n"
        << "  --threads <int>               Max worker threads (TBB/OMP). [default: auto]\n";
    printHelpOption();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        showUsage(argv[0]);
        return 1;
    }
    
    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    
    if (hasOption(opts, "--help") || filename.empty()) {
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }
    
    auto parallel = initParallelism(opts, false);
    initLogging("volt-coordination", parallel.threads);
    
    LammpsParser::Frame frame;
    if (!parseFrame(filename, frame)) return 1;
    
    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);
    
    CoordinationService analyzer;
    analyzer.setCutoff(getDouble(opts, "--cutoff", 3.2));
    analyzer.setRdfBins(getInt(opts, "--rdfBins", 500));
    
    spdlog::info("Starting coordination analysis...");
    json result = analyzer.compute(frame, outputBase);
    
    if (result.value("is_failed", false)) {
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }
    
    spdlog::info("Coordination analysis completed.");
    return 0;
}
