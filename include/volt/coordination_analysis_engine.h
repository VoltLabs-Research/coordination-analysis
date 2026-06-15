#pragma once

#include <volt/core/particle_property.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <map>
#include <string>

namespace Volt{

class CoordinationNumber{
public:
    CoordinationNumber(){}

    double cutoff() const{ return _cutoff; }
    void setCutoff(double newCutoff){ _cutoff = newCutoff; }

    const std::vector<double>& rdfX() const{ return _rdfX; }
    const std::vector<double>& rdY() const{ return _rdfY; }

    void transferComputationResults(class CoordinationAnalysisEngine* engine);

private:
    std::shared_ptr<ParticleProperty> _coordinationNumbers;
    double _cutoff;
    std::vector<double> _rdfX;
    std::vector<double> _rdfY;
};

// Per-type-pair RDF histogram entry
struct PartialRdfEntry {
    std::string pairType;   // e.g. "1-1", "1-2", "2-2"
    double binCenter;
    long long binCount;
};

class CoordinationAnalysisEngine{
public:
    CoordinationAnalysisEngine(
        ParticleProperty* positions,
        const SimulationCell& simCell,
        double cutoff,
        int rdfSampleCount,
        const std::vector<int>* types = nullptr
    )
        : _positions(positions),
          _simCell(simCell),
          _cutoff(cutoff),
          _rdfHistogram(rdfSampleCount, 0.0),
          _coordinationNumbers(new ParticleProperty(positions->size(), ParticleProperty::CoordinationProperty, 0, true)),
          _types(types){}

    void perform();

    ParticleProperty* positions() const{ return _positions; }
    const SimulationCell& cell() const{ return _simCell; }
    std::shared_ptr<ParticleProperty> coordinationNumbers() const{ return _coordinationNumbers; }
    double cutoff() const{ return _cutoff; }
    const std::vector<double>& rdfHistogram() const{ return _rdfHistogram; }
    const std::vector<PartialRdfEntry>& partialRdf() const{ return _partialRdf; }

    double _cutoff;
    SimulationCell _simCell;
    ParticleProperty* _positions;
    std::shared_ptr<ParticleProperty> _coordinationNumbers;
    std::vector<double> _rdfHistogram;

private:
    const std::vector<int>* _types;
    std::vector<PartialRdfEntry> _partialRdf;
};

}
