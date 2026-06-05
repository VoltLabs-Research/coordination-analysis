#include <volt/plugin/plugin_main.h>
#include <volt/plugin/analysis_plugin.h>
#include <volt/plugin/output_serializer.h>
#include <volt/coordination_analysis_engine.h>
#include <volt/core/frame_adapter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace Volt {

static std::array<double, 3> coordinationColor(int coordination, int minCoord, int maxCoord) {
    const double t = maxCoord > minCoord
        ? static_cast<double>(coordination - minCoord) / static_cast<double>(maxCoord - minCoord)
        : 0.5;
    return {
        0.15 + 0.75 * t,
        0.45 + 0.35 * (1.0 - std::abs(2.0 * t - 1.0)),
        0.90 - 0.70 * t
    };
}

class CoordinationPlugin : public Plugin::AnalysisPlugin<CoordinationPlugin> {
public:
    double cutoff = 3.2;
    int rdfBins = 500;

    std::string validate(const LammpsParser::Frame& frame) {
        if (!FrameAdapter::validateSimulationCell(frame.simulationCell))
            return "Invalid simulation cell";
        return "";
    }

    json run(const LammpsParser::Frame& frame,
             const std::shared_ptr<Particles::ParticleProperty>& positions,
             const std::string& outputBase) {
        CoordinationNumber coordNumber;
        coordNumber.setCutoff(cutoff);

        CoordinationAnalysisEngine engine(positions.get(), frame.simulationCell, cutoff, rdfBins);
        engine.perform();
        coordNumber.transferComputationResults(&engine);

        const auto& rdfX = coordNumber.rdfX();
        const auto& rdfY = coordNumber.rdY();
        auto coordProp = engine.coordinationNumbers();

        _coordinationValues.resize(static_cast<std::size_t>(frame.natoms), 0);
        _minCoordination = std::numeric_limits<int>::max();
        _maxCoordination = std::numeric_limits<int>::min();
        long long coordinationSum = 0;

        for (int i = 0; i < frame.natoms; ++i) {
            const int coord = coordProp ? coordProp->getInt(i) : 0;
            _coordinationValues[static_cast<std::size_t>(i)] = coord;
            _minCoordination = std::min(_minCoordination, coord);
            _maxCoordination = std::max(_maxCoordination, coord);
            coordinationSum += coord;
        }
        if (frame.natoms == 0) { _minCoordination = 0; _maxCoordination = 0; }

        json rdfRows = json::array();
        for (std::size_t i = 0; i < rdfX.size() && i < rdfY.size(); ++i) {
            rdfRows.push_back({{"r", rdfX[i]}, {"g_r", rdfY[i]}});
        }

        json result;
        result["main_listing"] = {
            {"cutoff", cutoff},
            {"rdf_bins", static_cast<int>(rdfRows.size())},
            {"total_atoms", frame.natoms},
            {"min_coordination", _minCoordination},
            {"max_coordination", _maxCoordination},
            {"mean_coordination", frame.natoms > 0
                ? static_cast<double>(coordinationSum) / static_cast<double>(frame.natoms) : 0.0}
        };
        result["sub_listings"] = {{"rdf", rdfRows}};

        if (!outputBase.empty()) {
            json chartWrapper;
            chartWrapper["export"]["ChartExporter"]["rdf"] = {
                {"r", rdfX}, {"g_r", rdfY}
            };
            const std::string chartPath = outputBase + "_rdf_chart.msgpack";
            if (JsonUtils::writeJsonMsgpackToFile(chartWrapper, chartPath, false))
                spdlog::info("RDF chart msgpack written to {}", chartPath);
        }

        return result;
    }

    void serialize(const LammpsParser::Frame& frame,
                   const std::shared_ptr<Particles::ParticleProperty>&,
                   const json& result,
                   const std::string& outputBase) {
        auto coordVals = _coordinationValues;
        int minC = _minCoordination, maxC = _maxCoordination;

        Plugin::serializePluginOutput(outputBase, frame, result, {
            .summaryFileSuffix = "_coordination",
            .bucketResolver = [&coordVals](std::size_t i) {
                return "Coordination_" + std::to_string(coordVals[i]);
            },
            .atomFieldWriter = [&coordVals, minC, maxC](MsgpackWriter& w, std::size_t i, int& count) {
                count = 3;
                const int coord = coordVals[i];
                const auto color = coordinationColor(coord, minC, maxC);
                w.write_key("coordination"); w.write_int(coord);
                w.write_key("coordination_color"); w.write_array_header(3);
                w.write_double(color[0]); w.write_double(color[1]); w.write_double(color[2]);
                w.write_key("color"); w.write_array_header(3);
                w.write_double(color[0]); w.write_double(color[1]); w.write_double(color[2]);
            }
        });
    }

private:
    std::vector<int> _coordinationValues;
    int _minCoordination = 0;
    int _maxCoordination = 0;
};

} // namespace Volt

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-coordination",
    .description = "Coordination Analysis",
    .options = {
        {"--cutoff", "float", "Cutoff radius for neighbor search", "3.2"},
        {"--rdf_bins", "int", "Number of bins for RDF calculation", "500"},
    }
};

VOLT_PLUGIN_MAIN(descriptor, [](const auto& opts, const Volt::LammpsParser::Frame& frame,
                                 const Volt::LammpsParser::Frame*, const std::string& outputBase) {
    Volt::CoordinationPlugin plugin;
    plugin.cutoff = Volt::CLI::getDouble(opts, "--cutoff", 3.2);
    plugin.rdfBins = Volt::CLI::getInt(opts, "--rdf_bins", 500);
    return plugin.compute(frame, outputBase);
})
