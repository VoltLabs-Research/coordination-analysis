#include <volt/coordination_analysis_engine.h>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <map>
#include <mutex>

namespace Volt{

void CoordinationAnalysisEngine::perform(){
    CutoffNeighborFinder neighborListBuilder;
    if(!neighborListBuilder.prepare(_cutoff, positions(), cell())){
        return;
    }

    const size_t particleCount = positions()->size();
    const double rdfBinSize = (_cutoff + EPSILON) / _rdfHistogram.size();
    int* coordOutput = _coordinationNumbers->dataInt();

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

    // Partial RDF per type-pair (serial, lightweight — bins are small maps)
    if(_types && _types->size() == particleCount){
        // Build unique type set
        std::map<int, bool> typeSet;
        for(int t : *_types) typeSet[t] = true;

        // Map from pair-key (ti*100000 + tj with ti<=tj) to histogram
        const int bins = static_cast<int>(_rdfHistogram.size());
        std::map<long long, std::vector<long long>> pairHist;
        for(auto& [ti, _a] : typeSet){
            for(auto& [tj, _b] : typeSet){
                if(ti <= tj){
                    long long key = (long long)ti * 100000LL + tj;
                    pairHist[key].assign(bins, 0LL);
                }
            }
        }

        for(size_t i = 0; i < particleCount; ++i){
            int ti = (*_types)[i];
            for(CutoffNeighborFinder::Query q(neighborListBuilder, i); !q.atEnd(); q.next()){
                int tj = (*_types)[q.current()];
                int ta = std::min(ti, tj);
                int tb = std::max(ti, tj);
                long long key = (long long)ta * 100000LL + tb;
                int bin = (int)(std::sqrt(q.distanceSquared()) / rdfBinSize);
                if(bin < bins) pairHist[key][bin]++;
            }
        }

        // Convert to PartialRdfEntry rows
        _partialRdf.clear();
        double stepSize = _cutoff / bins;
        for(auto& [key, hist] : pairHist){
            int ta = (int)(key / 100000LL);
            int tb = (int)(key % 100000LL);
            std::string pairName = std::to_string(ta) + "-" + std::to_string(tb);
            for(int b = 0; b < bins; ++b){
                _partialRdf.push_back({pairName, (b + 0.5) * stepSize, hist[b]});
            }
        }
    }
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
