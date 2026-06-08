#include <volt/coordination_analysis_engine.h>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Volt{

// Performs the actual computation. 
void CoordinationAnalysisEngine::perform(){
    // Prepare the neighbor list
    CutoffNeighborFinder neighborListBuilder;
    if(!neighborListBuilder.prepare(_cutoff, positions(), cell())){
        return;
    }

    const size_t particleCount = positions()->size();
    const double rdfBinSize = (_cutoff + EPSILON) / _rdfHistogram.size();
    int* coordOutput = _coordinationNumbers->dataInt();

    // Per-particle coordination is written at distinct indices (no hazard); the
    // RDF histogram is a reduction: each range accumulates a thread-local
    // histogram and the join sums them element-wise. Replaces the previous
    // manual std::thread chunking + mutex (which did not compose with the TBB
    // thread pool governed by --threads).
    _rdfHistogram = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, particleCount),
        std::vector<double>(_rdfHistogram.size(), 0.0),
        [&](const tbb::blocked_range<size_t>& r, std::vector<double> localRDF){
            for(size_t i = r.begin(); i < r.end(); ++i){
                int coordNumber = 0;
                for(CutoffNeighborFinder::Query q(neighborListBuilder, i); !q.atEnd(); q.next()){
                    coordNumber++;
                    size_t bin = (size_t)(std::sqrt(q.distanceSquared()) / rdfBinSize);
                    if(bin < localRDF.size()) localRDF[bin]++;
                }
                coordOutput[i] = coordNumber;
            }
            return localRDF;
        },
        [](std::vector<double> a, const std::vector<double>& b){
            for(size_t i = 0; i < a.size(); ++i) a[i] += b[i];
            return a;
        });
}

void CoordinationNumber::transferComputationResults(CoordinationAnalysisEngine* engine){
    CoordinationAnalysisEngine* eng = static_cast<CoordinationAnalysisEngine*>(engine);
    _coordinationNumbers = eng->coordinationNumbers();

    const auto& hist = eng->rdfHistogram();

    _rdfY.resize(hist.size());
    _rdfX.resize(hist.size());

    double rho = static_cast<double>(eng->positions()->size()) / eng->cell().volume3D();
    double constant = 4.0 / 3.0 * PI * rho * eng->positions()->size();
    double stepSize = eng->cutoff() / _rdfX.size();

    for(size_t i = 0; i < _rdfX.size(); i++){
        double r  = stepSize * i;
        double r2 = r + stepSize;
        _rdfX[i] = r + 0.5 * stepSize;
        _rdfY[i] = hist[i] / (constant * (r2 * r2 * r2 - r * r * r));
    }
}
}
