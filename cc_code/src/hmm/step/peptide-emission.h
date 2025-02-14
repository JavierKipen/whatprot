/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_PEPTIDE_EMISSION_H
#define WHATPROT_HMM_STEP_PEPTIDE_EMISSION_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "common/radiometry.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/peptide-step.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/kd-range.h"

namespace whatprot {

class PeptideEmission : public PeptideStep {
public:
    PeptideEmission(const Radiometry& radiometry,
                    unsigned int timestep,
                    int max_num_dyes,
                    const SequencingModel& seq_model,
                    const SequencingSettings& seq_settings);
    PeptideEmission(const PeptideEmission& other);
    virtual ~PeptideEmission();
    double& prob(int channel, int num_dyes);
    double prob(int channel, int num_dyes) const;
    virtual void prune_forward(KDRange* range, bool* allow_detached) override;
    virtual void prune_backward(KDRange* range, bool* allow_detached) override;
    PeptideStateVector* forward_or_backward(const PeptideStateVector& input,
                                            unsigned int* num_edmans) const;
    virtual PeptideStateVector* forward(
            const PeptideStateVector& input,
            unsigned int* num_edmans) const override;
    virtual PeptideStateVector* backward(
            const PeptideStateVector& input,
            unsigned int* num_edmans) const override;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
    const Radiometry& radiometry;
    unsigned int timestep;
    KDRange pruned_range;
    bool allow_detached;
    std::vector<double>* values;
    bool i_am_a_copy;
    unsigned int num_channels;
    int max_num_dyes;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_PEPTIDE_EMISSION_H
