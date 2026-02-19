#include <volt/coordination_service.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <spdlog/spdlog.h>

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
    std::vector<int> coord(frame.natoms);
    for(int i = 0; i < frame.natoms; i++){
        coord[i] = coordProp->getInt(i);
    }

    json result = AnalysisResult::success();
    AnalysisResult::addTiming(result, startTime);
    result["cutoff"] = _cutoff;
    result["rdf"]["x"] = rdfX;
    result["rdf"]["y"] = rdfY;
    result["coordination"] = coord;

    if(!outputFile.empty()){
        // TODO: Implement msgpack export in standalone package
        spdlog::warn("File output not yet implemented in standalone package");
    }

    return result;
}

}
