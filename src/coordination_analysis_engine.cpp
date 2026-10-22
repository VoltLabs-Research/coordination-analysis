#include <volt/coordination_analysis_engine.h>
#include <volt/core/progress.h>
#include <spdlog/spdlog.h>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <map>
#include <mutex>

namespace Volt{

void CoordinationAnalysisEngine::perform(){
    spdlog::info("Building neighbor lists (cutoff {:.5g}) for {} atoms...", _cutoff, positions()->size());
    Stopwatch timer;
    CutoffNeighborFinder neighborListBuilder;
    if(!neighborListBuilder.prepare(_cutoff, positions(), cell())){
        spdlog::error("Failed to build neighbor lists; coordination numbers and RDF are left at zero");
        return;
    }
    spdlog::info("Neighbor lists ready in {}", timer.elapsed());

    const size_t particleCount = positions()->size();
    const double rdfBinSize = (_cutoff + EPSILON) / _rdfHistogram.size();
    int* coordOutput = _coordinationNumbers->dataInt();

    spdlog::info("Computing coordination numbers and RDF for {} atoms...", particleCount);
    ProgressReporter progress("Computing coordination numbers", particleCount);
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
            progress.advance(r.size());
            return localRDF;
        },
        [](std::vector<double> a, const std::vector<double>& b){
            for(size_t i = 0; i < a.size(); ++i) a[i] += b[i];
            return a;
        });
    spdlog::info("Computed coordination numbers and RDF for {} atoms in {}", particleCount, progress.elapsed());

    if(_types && _types->size() == particleCount){
        timer.restart();
        std::map<int, bool> typeSet;
        for(int t : *_types) typeSet[t] = true;

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

        spdlog::info("Computing partial RDFs for {} atom types ({} pairs)...", typeSet.size(), pairHist.size());
        ProgressReporter pairProgress("Computing partial RDFs", particleCount);
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
            pairProgress.update(i + 1);
        }

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
        spdlog::info("Computed partial RDFs in {}", timer.elapsed());
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
