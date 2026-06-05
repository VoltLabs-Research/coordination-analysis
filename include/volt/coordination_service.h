#pragma once

#include <volt/coordination_analysis_engine.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt {

using json = nlohmann::json;

class CoordinationService {
public:
    void setCutoff(double v);
    void setRdfBins(int v);

    json compute(const LammpsParser::Frame& frame, const std::string& outputBase);

private:
    double _cutoff = 3.2;
    int _rdfBins = 500;
};

} // namespace Volt
