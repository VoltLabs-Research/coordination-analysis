#include <volt/plugin/plugin_entry.h>
#include <volt/coordination_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = CoordinationService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--cutoff", "Cutoff radius for neighbor search", 3.2, &S::setCutoff),
    opt("--rdf_bins", "Number of bins for RDF calculation", 500, &S::setRdfBins),
};

VOLT_SERVICE_PLUGIN("volt-coordination", "Coordination Analysis", S, bindings)
