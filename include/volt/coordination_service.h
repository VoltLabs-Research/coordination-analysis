#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <volt/core/particle_property.h>
#include <volt/coordination_analysis_engine.h>
#include <string>

namespace Volt{
using json = nlohmann::json;

class CoordinationService{
public:
    CoordinationService();

    void setCutoff(double cutoff);
    void setRdfBins(int bins);

    json compute(
        const LammpsParser::Frame &frame,
        const std::string &outputFilename = ""
    );

private:
    double _cutoff;
    int _rdfBins;
};

}
