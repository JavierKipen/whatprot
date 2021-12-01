/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "binomial-transition.h"

// Standard C++ library headers:
#include <algorithm>
#include <limits>

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/parameter-fitter.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace whatprot {

namespace {
using std::numeric_limits;
}  // namespace

BinomialTransition::BinomialTransition(double q, int channel)
        : q(q), channel(channel) {
    length = 1;
    size = 1;
    values.resize(size);
    values[0] = 1.0;
}

void BinomialTransition::reserve(unsigned int max_n) {
    if (max_n + 1 <= length) {
        return;
    }
    unsigned int prev_length = length;
    length = max_n + 1;
    size = length * (length + 1) / 2;
    values.resize(size);
    double p = (double)1 - q;
    for (unsigned int i = prev_length; i < length; i++) {
        prob(i, 0) = prob(i - 1, 0) * q;
        for (unsigned int j = 1; j < i; j++) {
            prob(i, j) = prob(i - 1, j) * q + prob(i - 1, j - 1) * p;
        }
        prob(i, i) = prob(i - 1, i - 1) * p;
    }
}

double& BinomialTransition::prob(unsigned int from, unsigned int to) {
    return values[from * (from + 1) / 2 + to];
}

double BinomialTransition::prob(unsigned int from, unsigned int to) const {
    return values[from * (from + 1) / 2 + to];
}

void BinomialTransition::prune_forward(KDRange* range, bool* allow_detached) {
    forward_range = *range;
    range->min[1 + channel] = 0;
    backward_range = *range;
}

void BinomialTransition::prune_backward(KDRange* range, bool* allow_detached) {
    backward_range = backward_range.intersect(*range);
    *range = backward_range;
    range->max[1 + channel] = numeric_limits<unsigned int>::max();
    forward_range = forward_range.intersect(*range);
    *range = forward_range;
}

void BinomialTransition::forward(unsigned int* num_edmans,
                                 PeptideStateVector* psv) const {
    TensorVectorIterator* itr =
            psv->tensor.vector_iterator(forward_range, 1 + channel);
    while (!itr->done()) {
        this->forward(itr->get());
        itr->advance();
    }
    psv->range = backward_range;
}

void BinomialTransition::forward(Vector* v) const {
    unsigned int to_min = backward_range.min[1 + channel];
    unsigned int to_max = backward_range.max[1 + channel];
    for (unsigned int to = to_min; to < to_max; to++) {
        double v_to = 0.0;
        unsigned int from_min = std::max(to, forward_range.min[1 + channel]);
        unsigned int from_max = forward_range.max[1 + channel];
        for (unsigned int from = from_min; from < from_max; from++) {
            v_to += prob(from, to) * (*v)[from];
        }
        (*v)[to] = v_to;
    }
}

void BinomialTransition::backward(const PeptideStateVector& input,
                                  unsigned int* num_edmans,
                                  PeptideStateVector* output) const {
    ConstTensorVectorIterator* in_itr =
            input.tensor.const_vector_iterator(backward_range, 1 + channel);
    TensorVectorIterator* out_itr =
            output->tensor.vector_iterator(forward_range, 1 + channel);
    while (!out_itr->done()) {
        this->backward(*in_itr->get(), out_itr->get());
        in_itr->advance();
        out_itr->advance();
    }
    output->range = forward_range;
}

void BinomialTransition::backward(const Vector& input, Vector* output) const {
    int from_min = forward_range.min[1 + channel];
    int from_max = forward_range.max[1 + channel];
    for (int from = from_max - 1; from >= from_min; from--) {
        double v_from = 0.0;
        int to_min = backward_range.min[1 + channel];
        int to_max = std::min(from, (int)backward_range.max[1 + channel] - 1);
        for (int to = to_min; to <= to_max; to++) {
            v_from += prob(from, to) * input[to];
        }
        (*output)[from] = v_from;
    }
}

void BinomialTransition::improve_fit(
        const PeptideStateVector& forward_psv,
        const PeptideStateVector& backward_psv,
        const PeptideStateVector& next_backward_psv,
        unsigned int num_edmans,
        double probability,
        ParameterFitter* fitter) const {
    ConstTensorVectorIterator* f_itr = forward_psv.tensor.const_vector_iterator(
            forward_range, 1 + channel);
    ConstTensorVectorIterator* b_itr =
            backward_psv.tensor.const_vector_iterator(forward_range,
                                                      1 + channel);
    ConstTensorVectorIterator* nb_itr =
            next_backward_psv.tensor.const_vector_iterator(backward_range,
                                                           1 + channel);
    while (!f_itr->done()) {
        this->improve_fit(*f_itr->get(),
                          *b_itr->get(),
                          *nb_itr->get(),
                          probability,
                          fitter);
        f_itr->advance();
        b_itr->advance();
        nb_itr->advance();
    }
}

void BinomialTransition::improve_fit(const Vector& forward_vector,
                                     const Vector& backward_vector,
                                     const Vector& next_backward_vector,
                                     double probability,
                                     ParameterFitter* fitter) const {
    int from_min = forward_range.min[1 + channel];
    int from_max = forward_range.max[1 + channel];
    // Note that we can ignore when starting location (from) is 0 because then
    // there are no dyes, it's irrelevant. We would be adding 0s.
    for (int from = from_max - 1; from >= std::max(from_min, 1); from--) {
        double p_state =
                forward_vector[from] * backward_vector[from] / probability;
        fitter->denominator += p_state * (double)from;
        int to_min = backward_range.min[1 + channel];
        int to_max = backward_range.max[1 + channel];
        // We can ignore when to and from are equal, because no dyes are lost
        // then, so it gives us nothing else for the numerator; we would be
        // adding zero.
        for (int to = to_min; to < std::min(to_max, from); to++) {
            double p_transition = forward_vector[from] * prob(from, to)
                                  * next_backward_vector[to] / probability;
            fitter->numerator += p_transition * (double)(from - to);
        }
    }
}

}  // namespace whatprot
