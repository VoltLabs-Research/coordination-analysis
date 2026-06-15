#include <volt/coordination_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/plugin/output_serializer.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace Volt {

void CoordinationService::setCutoff(double v) { _cutoff = v; }
void CoordinationService::setRdfBins(int v) { _rdfBins = v; }

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

json CoordinationService::compute(const LammpsParser::Frame& frame, const std::string& outputBase) {
    if (frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");
    if (!FrameAdapter::validateSimulationCell(frame.simulationCell))
        return AnalysisResult::failure("Invalid simulation cell");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if (!positions) return AnalysisResult::failure("Failed to create position property");

    CoordinationNumber coordNumber;
    coordNumber.setCutoff(_cutoff);

    const std::vector<int>* typesPtr = frame.types.empty() ? nullptr : &frame.types;
    CoordinationAnalysisEngine engine(positions.get(), frame.simulationCell, _cutoff, _rdfBins, typesPtr);
    engine.perform();
    coordNumber.transferComputationResults(&engine);

    const auto& rdfX = coordNumber.rdfX();
    const auto& rdfY = coordNumber.rdY();
    auto coordProp = engine.coordinationNumbers();

    std::vector<int> coordinationValues(static_cast<std::size_t>(frame.natoms), 0);
    int minCoordination = std::numeric_limits<int>::max();
    int maxCoordination = std::numeric_limits<int>::min();
    long long coordinationSum = 0;

    for (int i = 0; i < frame.natoms; ++i) {
        const int coord = coordProp ? coordProp->getInt(i) : 0;
        coordinationValues[static_cast<std::size_t>(i)] = coord;
        minCoordination = std::min(minCoordination, coord);
        maxCoordination = std::max(maxCoordination, coord);
        coordinationSum += coord;
    }
    if (frame.natoms == 0) { minCoordination = 0; maxCoordination = 0; }

    json rdfRows = json::array();
    for (std::size_t i = 0; i < rdfX.size() && i < rdfY.size(); ++i) {
        rdfRows.push_back({{"r", rdfX[i]}, {"g_r", rdfY[i]}});
    }

    json result;
    result["main_listing"] = {
        {"cutoff", _cutoff},
        {"rdf_bins", static_cast<int>(rdfRows.size())},
        {"total_atoms", frame.natoms},
        {"min_coordination", minCoordination},
        {"max_coordination", maxCoordination},
        {"mean_coordination", frame.natoms > 0
            ? static_cast<double>(coordinationSum) / static_cast<double>(frame.natoms) : 0.0}
    };
    result["sub_listings"] = {{"rdf", rdfRows}};

    if (!outputBase.empty()) {
        json chartWrapper;
        chartWrapper["export"]["ChartExporter"]["rdf"] = {
            {"r", rdfX}, {"g_r", rdfY}
        };
        const std::string chartPath = outputBase + "_rdf_chart.parquet";
        if (JsonUtils::writeJsonToParquet(chartWrapper, chartPath))
            spdlog::info("RDF chart parquet written to {}", chartPath);

        // Partial RDF histogram (per type-pair)
        const auto& partialRdf = engine.partialRdf();
        if(!partialRdf.empty()){
            json histWrapper;
            auto& rows = histWrapper["export"]["data"] = json::array();
            for(const auto& e : partialRdf){
                rows.push_back({{"pair_type", e.pairType}, {"bin_center", e.binCenter}, {"bin_count", e.binCount}});
            }
            const std::string histPath = outputBase + "_rdf_histogram.parquet";
            if(JsonUtils::writeJsonToParquet(histWrapper, histPath))
                spdlog::info("Partial RDF histogram parquet written to {}", histPath);
        }

        int minC = minCoordination, maxC = maxCoordination;
        Plugin::serializePluginOutput(outputBase, frame, result, {
            .summaryFileSuffix = "_coordination",
            .bucketResolver = [&coordinationValues](std::size_t i) {
                return "Coordination_" + std::to_string(coordinationValues[i]);
            },
            .perAtomColumnWriter = [&coordinationValues, minC, maxC](ColumnarAtomWriter& w, std::size_t i) {
                const int coord = coordinationValues[i];
                const auto color = coordinationColor(coord, minC, maxC);
                w.field("coordination", coord);
                w.field("coordination_color", std::vector<double>{color[0], color[1], color[2]});
                w.field("color", std::vector<double>{color[0], color[1], color[2]});
            }
        });
    }

    return result;
}

} // namespace Volt
