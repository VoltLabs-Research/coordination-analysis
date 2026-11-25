#pragma once

#include <volt/coordination_analysis_engine.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <array>
#include <string>

namespace Volt {

using json = nlohmann::json;

class CoordinationService {
public:
    void setCutoff(double v);
    void setRdfBins(int v);

    json compute(const LammpsParser::Frame& frame, const std::string& outputBase);

private:
    // The color of an atom in the coordination range of the frame: blue for the lowest, red for the highest.
    static std::array<double, 3> coordinationColor(int coordination, int minCoord, int maxCoord);

    double _cutoff = 3.2;
    int _rdfBins = 500;
};

}
