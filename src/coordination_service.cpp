#include <volt/coordination_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

CoordinationService::CoordinationService()
    : _cutoff(3.2),
      _rdfBins(500){}

void CoordinationService::setCutoff(double cutoff){
    _cutoff = cutoff;
}

void CoordinationService::setRdfBins(int bins){
    _rdfBins = bins;
}

namespace{

std::array<double, 3> coordinationColor(int coordination, int minCoordination, int maxCoordination){
    const double t = maxCoordination > minCoordination
        ? static_cast<double>(coordination - minCoordination) / static_cast<double>(maxCoordination - minCoordination)
        : 0.5;
    return {
        0.15 + 0.75 * t,
        0.45 + 0.35 * (1.0 - std::abs(2.0 * t - 1.0)),
        0.90 - 0.70 * t
    };
}

}

json CoordinationService::compute(const LammpsParser::Frame &frame, const std::string& outputFile){
    auto startTime = std::chrono::high_resolution_clock::now();

    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    if(!FrameAdapter::validateSimulationCell(frame.simulationCell))
        return AnalysisResult::failure("Invalid simulation cell");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    spdlog::info("Starting coordination analysis (cutoff = {}, bins = {})...", _cutoff, _rdfBins);
    CoordinationNumber coordNumber;
    coordNumber.setCutoff(_cutoff);

    CoordinationAnalysisEngine engine(
        positions.get(),
        frame.simulationCell,
        _cutoff,
        _rdfBins
    );

    engine.perform(),
    coordNumber.transferComputationResults(&engine);

    const auto &rdfX = coordNumber.rdfX();
    const auto &rdfY = coordNumber.rdY();

    auto coordProp = engine.coordinationNumbers();

    std::vector<int> coordinationValues(static_cast<std::size_t>(frame.natoms), 0);
    int minCoordination = std::numeric_limits<int>::max();
    int maxCoordination = std::numeric_limits<int>::min();
    long long coordinationSum = 0;
    for(int i = 0; i < frame.natoms; i++){
        const int coord = coordProp ? coordProp->getInt(i) : 0;
        coordinationValues[static_cast<std::size_t>(i)] = coord;
        minCoordination = std::min(minCoordination, coord);
        maxCoordination = std::max(maxCoordination, coord);
        coordinationSum += coord;
    }
    if(frame.natoms == 0){
        minCoordination = 0;
        maxCoordination = 0;
    }

    json rdfRows = json::array();
    for(std::size_t i = 0; i < rdfX.size() && i < rdfY.size(); ++i){
        rdfRows.push_back({
            {"r", rdfX[i]},
            {"g_r", rdfY[i]}
        });
    }

    json result;
    result["main_listing"] = {
        { "cutoff", _cutoff },
        { "rdf_bins", static_cast<int>(rdfRows.size()) },
        { "total_atoms", frame.natoms },
        { "min_coordination", minCoordination },
        { "max_coordination", maxCoordination },
        { "mean_coordination", frame.natoms > 0 ? static_cast<double>(coordinationSum) / static_cast<double>(frame.natoms) : 0.0 }
    };
    result["sub_listings"] = {
        { "rdf", rdfRows }
    };

    json perAtom = json::array();
    for(int i = 0; i < frame.natoms; i++){
        const int coord = coordinationValues[static_cast<std::size_t>(i)];
        const auto color = coordinationColor(coord, minCoordination, maxCoordination);
        perAtom.push_back({
            { "id", frame.ids[i] },
            { "coordination", coord },
            { "coordination_color", {color[0], color[1], color[2]} }
        });
    }
    result["per-atom-properties"] = perAtom;

    if(!outputFile.empty()){
        const std::string outputPath = outputFile + "_coordination.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Coordination msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write coordination msgpack: {}", outputPath);
        }

        // --- atoms.msgpack (AtomisticExporter) ---
        // Canonical per-atom envelope with the coordination number exposed
        // as an extra column (OVITO's CoordinationAnalysisModifier
        // publishes a `CoordinationProperty` scalar per atom). Grouping
        // atoms by coordination bucket gives the viewport a natural
        // categorical axis for colouring.
        json atomsByCoord;
        json atomProperties = json::array();
        std::map<int, int> coordCounts;
        for(int i = 0; i < frame.natoms; i++){
            const Point3& pos = frame.positions[i];
            const int coord = coordinationValues[static_cast<std::size_t>(i)];
            const auto color = coordinationColor(coord, minCoordination, maxCoordination);
            const std::string bucket = "Coordination_" + std::to_string(coord);
            const int atomId = i < static_cast<int>(frame.ids.size()) ? frame.ids[i] : i;
            coordCounts[coord]++;
            json atom = {
                {"id", atomId},
                {"pos", {pos.x(), pos.y(), pos.z()}},
                {"structure_id", coord},
                {"structure_name", bucket},
                {"cluster_id", 0},
                {"coordination", coord},
                {"coordination_color", {color[0], color[1], color[2]}},
                {"color", {color[0], color[1], color[2]}}
            };
            atomsByCoord[bucket].push_back(atom);
            atomProperties.push_back({
                {"id", atomId},
                {"coordination", coord},
                {"structure_name", bucket}
            });
        }
        json structuresListing = json::array();
        for(const auto& [coord, count] : coordCounts){
            structuresListing.push_back({
                {"structure_id", coord},
                {"structure_name", "Coordination_" + std::to_string(coord)},
                {"atom_count", count}
            });
        }
        json exportWrapper;
        exportWrapper["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"structure_count", static_cast<int>(coordCounts.size())}
        };
        exportWrapper["sub_listings"] = { {"structures", structuresListing} };
        exportWrapper["per-atom-properties"] = atomProperties;
        exportWrapper["export"] = json::object();
        exportWrapper["export"]["AtomisticExporter"] = atomsByCoord;
        const std::string atomsPath = outputFile + "_atoms.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(exportWrapper, atomsPath, false)){
            spdlog::info("Exported atoms data to: {}", atomsPath);
        }else{
            spdlog::warn("Could not write atoms msgpack: {}", atomsPath);
        }

        json chartWrapper;
        chartWrapper["export"]["ChartExporter"]["rdf"] = {
            {"r", rdfX},
            {"g_r", rdfY}
        };
        const std::string chartPath = outputFile + "_rdf_chart.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(chartWrapper, chartPath, false)){
            spdlog::info("RDF chart msgpack written to {}", chartPath);
        }else{
            spdlog::warn("Could not write RDF chart msgpack: {}", chartPath);
        }
    }

    return result;
}

}
