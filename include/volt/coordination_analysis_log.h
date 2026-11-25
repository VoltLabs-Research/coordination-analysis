#pragma once

#include <volt/core/run_log.h>
#include <volt/coordination_analysis_engine.h>

#include <vector>

namespace Volt{

// The execution log of one analysis: each phase hands it what it built.
class CoordinationAnalysisLog : public RunLog{
public:
    using RunLog::dump;

    // The engine fills them together with the RDF histogram, so their line covers both.
    void dump(const ParticleProperty& coordinationNumbers);
    void dump(const std::vector<PartialRdfEntry>& partialRdf);
    // The coordination numbers and the normalized RDF, once the engine handed them over.
    void dump(const CoordinationNumber& coordination);
};

}
