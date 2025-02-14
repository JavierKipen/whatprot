/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "kd-range.h"

// Standard C++ library headers:
#include <algorithm>
#include <vector>

namespace whatprot {

KDRange KDRange::intersect(const KDRange& other) const {
    KDRange result;
    result.min.resize(min.size());
    result.max.resize(max.size());
    for (unsigned int i = 0; i < result.min.size(); i++) {
        result.min[i] = std::max(min[i], other.min[i]);
        result.max[i] = std::min(max[i], other.max[i]);
    }
    return result;
}

bool KDRange::is_empty() const {
    for (unsigned int i = 0; i < min.size(); i++) {
        if (min[i] >= max[i]) {
            return true;
        }
    }
    return false;
}

bool KDRange::includes_zero() const {
    for (unsigned int i = 0; i < min.size(); i++) {
        if (min[i] > 0) {
            return false;
        }
        if (max[i] <= 0) {
            return false;
        }
    }
    return true;
}

}  // namespace whatprot
