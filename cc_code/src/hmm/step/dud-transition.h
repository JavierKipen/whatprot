/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_DUD_TRANSITION_H
#define WHATPROT_HMM_STEP_DUD_TRANSITION_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/binomial-transition.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

class DudTransition : public BinomialTransition {
public:
    DudTransition(double q, int channel);
    DudTransition(const DudTransition& other);
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_DUD_TRANSITION_H