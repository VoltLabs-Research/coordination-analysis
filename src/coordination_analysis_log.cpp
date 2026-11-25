#include <volt/coordination_analysis_log.h>
#include <volt/core/progress.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstdint>
#include <limits>

namespace Volt{

void CoordinationAnalysisLog::dump(const ParticleProperty& coordinationNumbers){
    spdlog::info("Computed coordination numbers and RDF for {} atoms in {}", coordinationNumbers.size(), phaseElapsed());
}

void CoordinationAnalysisLog::dump(const std::vector<PartialRdfEntry>&){
    spdlog::info("Computed partial RDFs in {}", phaseElapsed());
}

void CoordinationAnalysisLog::dump(const CoordinationNumber& coordination){
    const ParticleProperty& numbers = *coordination.coordinationNumbers();
    const std::size_t atomCount = numbers.size();
    int minCoordination = std::numeric_limits<int>::max();
    int maxCoordination = std::numeric_limits<int>::min();
    long long coordinationSum = 0;
    for(std::size_t i = 0; i < atomCount; ++i){
        const int coord = numbers.getInt(i);
        minCoordination = std::min(minCoordination, coord);
        maxCoordination = std::max(maxCoordination, coord);
        coordinationSum += coord;
    }

    std::vector<std::size_t> coordinationHistogram(static_cast<std::size_t>(maxCoordination - minCoordination) + 1, 0);
    for(std::size_t i = 0; i < atomCount; ++i)
        ++coordinationHistogram[static_cast<std::size_t>(numbers.getInt(i) - minCoordination)];
    const auto mostCommon = std::max_element(coordinationHistogram.begin(), coordinationHistogram.end());
    spdlog::info("Coordination numbers: min {}, mean {:.2f}, max {}, most common {} ({:.1f}% of atoms)",
        minCoordination, static_cast<double>(coordinationSum) / static_cast<double>(atomCount), maxCoordination,
        minCoordination + static_cast<int>(mostCommon - coordinationHistogram.begin()),
        percentOf(*mostCommon, static_cast<std::uint64_t>(atomCount)));

    const auto& rdfX = coordination.rdfX();
    const auto& rdfY = coordination.rdY();
    if(!rdfY.empty() && rdfX.size() == rdfY.size()){
        const auto peak = std::max_element(rdfY.begin(), rdfY.end());
        spdlog::info("RDF: {} bins up to r = {:.5g}, peak g(r) = {:.3f} at r = {:.3f}",
            rdfY.size(), coordination.cutoff(), *peak, rdfX[static_cast<std::size_t>(peak - rdfY.begin())]);
    }
}

}
