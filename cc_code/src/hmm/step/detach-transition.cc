/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "detach-transition.h"

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

DetachTransition::DetachTransition(double p_detach, unsigned int num_channels)
        : p_detach(p_detach), num_channels(num_channels), omit_detach(false) {}

void DetachTransition::prune_forward(KDRange* range, bool* allow_detached) {
    forward_range = *range;
    *allow_detached = true;
}

void DetachTransition::prune_backward(KDRange* range, bool* allow_detached) {
    // Pruned backwards range is always the intersection of the provided forward
    // and backwards ranges, regardless of whether 0 is in that range. Note,
    // this may not hold if forward range is not guaranteed to include 0 - which
    // will happen if DetachTransition is not used directly after
    // Bleach Transition objects covering all channels.
    backward_range = forward_range.intersect(*range);
    // If a detached peptide is possible, then the pruned forward range should
    // remain unchanged from what was provided, since any value in the forward
    // range can transition to a detached peptide. Otherwise the forward and
    // backward ranges should be identical.
    if (!*allow_detached) {
        forward_range = backward_range;
        omit_detach = true;
    }
    // No matter what how it was computed, we need to communicate what the
    // forward range ended up being.
    *range = forward_range;
}

void DetachTransition::forward(unsigned int* num_edmans,
                               PeptideStateVector* psv) const {
    double sum = 0.0;
    // Note, forward range is always strictly larger than backwards range.
    TensorIterator* itr = psv->tensor.iterator(forward_range);
    while (!itr->done()) {
        double value = *itr->get();
        // We will go ahead and modify this even if it is not in the output
        // range, because (1) if it's not read, it doesn't matter what's there,
        // and (2) this is likely far faster (though it hasn't been tested) than
        // having a messy if statement conditioned on whether this location is
        // in the backwards range.
        *itr->get() = value * (1 - p_detach);
        sum += value;
        itr->advance();
    }
    psv->tensor.values[(*num_edmans) * psv->tensor.strides[0]] +=
            p_detach * sum;
    psv->range = backward_range;
}

void DetachTransition::backward(const PeptideStateVector& input,
                                unsigned int* num_edmans,
                                PeptideStateVector* output) const {
    // If there is no detachment event in the (backward) input, this is much
    // simpler, since detachment from literally any position need not be
    // considered the way it is when detachment is possible. This case resembles
    // the "Second pass" section below, with one key difference; an equals sign
    // instead of a plus equals sign.
    if (omit_detach) {
        ConstTensorIterator* in_itr = input.tensor.const_iterator(backward_range);
        TensorIterator* out_itr_2 = output->tensor.iterator(backward_range);
        while (!out_itr_2->done()) {
            // Here we use an equals sign, NOT the same as the similar code
            // below under "Second pass."
            *out_itr_2->get() = (1 - p_detach) * (*in_itr->get());
            in_itr->advance();
            out_itr_2->advance();
        }
    } else {
        // Easiest way to handle potentially mismatched forward and backward
        // ranges is to do this in two passes. Note that we can't just ignore
        // this like we do for the forward() function, because we don't want to
        // accidentally process garbage input.
        //
        // First pass.
        TensorIterator* out_itr_1 = output->tensor.iterator(forward_range);
        while (!out_itr_1->done()) {
            *out_itr_1->get() = p_detach * input.p_detached;
            out_itr_1->advance();
        }
        // Second pass.
        ConstTensorIterator* in_itr = input.tensor.const_iterator(backward_range);
        TensorIterator* out_itr_2 = output->tensor.iterator(backward_range);
        while (!out_itr_2->done()) {
            // Here we use a plus equals sign, NOT the same as the similar code
            // above under "if (omit_detach)".
            *out_itr_2->get() += (1 - p_detach) * (*in_itr->get());
            in_itr->advance();
            out_itr_2->advance();
        }
    }
    // And finish it up....
    output->range = forward_range;
}

void DetachTransition::improve_fit(const PeptideStateVector& forward_psv,
                                   const PeptideStateVector& backward_psv,
                                   const PeptideStateVector& next_backward_psv,
                                   unsigned int num_edmans,
                                   double probability,
                                   SequencingModelFitter* fitter) const {
    double forward_sum = 0.0;
    double forward_backward_sum = 0.0;
    ConstTensorIterator* f_itr = forward_psv.tensor.const_iterator(
        forward_range);
    ConstTensorIterator* b_itr = backward_psv.tensor.const_iterator(
        forward_range);
    int t_idx = -1;  // used to skip 0th entry of every timestep (see below)
    while(!f_itr->done()) {
        // Skip the zeroth value in every timestep because this is the entry for
        // zero of every dye color. These don't provide evidence of detachment
        // one way or the other. However, we only have to do this if the range
        // includes the detached state; omit_detach tells us this.
        if (!omit_detach && t_idx != f_itr->loc[0]) {
            t_idx = f_itr->loc[0];
            f_itr->advance();
            b_itr->advance();
            continue;
        }
        forward_sum += *f_itr->get();
        forward_backward_sum += (*f_itr->get()) * (*b_itr->get());
        f_itr->advance();
        b_itr->advance();
    }
    // And continue....
    if (!omit_detach) {
        fitter->p_detach_fit.numerator +=
                forward_sum * p_detach
                * next_backward_psv.p_detached
                / probability;
    }
    // Probability of being in a state that can detach is 1.0, because all
    // states can detach (although we are ignoring the case where there are no
    // amino acids left but this probably shouldn't cause any serious issues).
    fitter->p_detach_fit.denominator += forward_backward_sum / probability;
}

}  // namespace whatprot
