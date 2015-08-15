#ifndef che_molecule_h
#define che_molecule_h

#include <boost/graph/undirected_graph.hpp>
#include "che/atom.h"
#include "che/bond.h"

namespace biosim {
  namespace che {
    // molecule is a graph
    using molecule = boost::undirected_graph<che::atom, che::bond>;
  } // namespace che
} // namespace biosim

#endif // che_molecule_h
