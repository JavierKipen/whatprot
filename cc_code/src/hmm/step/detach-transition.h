/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_DETACH_TRANSITION_H
#define WHATPROT_HMM_STEP_DETACH_TRANSITION_H

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/peptide-step.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "util/kd-range.h"

namespace whatprot {

class DetachTransition : public PeptideStep {
public:
    DetachTransition(double p_detach);
    virtual void prune_forward(KDRange* range) override;
    virtual void prune_backward(KDRange* range) override;
    virtual void forward(unsigned int* num_edmans,
                         PeptideStateVector* psv) const override;
    // The 'input' and 'output' parameters MUST NOT be the same
    // PeptideStateVector, if they are, this will not do what you want.
    virtual void backward(const PeptideStateVector& input,
                          unsigned int* num_edmans,
                          PeptideStateVector* output) const override;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
    KDRange forward_range;
    KDRange backward_range;
    double p_detach;
    unsigned int num_channels;
    bool omit_detach;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_DETACH_TRANSITION_H
