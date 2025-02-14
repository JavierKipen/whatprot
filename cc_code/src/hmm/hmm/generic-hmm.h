/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_HMM_GENERIC_HMM_H
#define WHATPROT_HMM_HMM_GENERIC_HMM_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/step/step.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

// V is the state vector type.
// S is the subclass of the Step class being used.
template <typename V, typename S>
class GenericHMM {
public:
    GenericHMM(unsigned int num_timesteps) : num_timesteps(num_timesteps) {}

    ~GenericHMM() {
        for (S* s : steps) {
            delete s;
        }
    }

    virtual V* create_states_forward() const = 0;

    virtual V* create_states_backward() const = 0;

    // This computes the probability of the provided dye seq producing the
    // provided radiometry. To do this efficiently, it uses a modified version
    // of the forward algorithm.
    virtual double probability() const {
        unsigned int num_edmans = 0;
        auto step = steps.begin();  // const_iterator type
        V* states_in = create_states_forward();
        states_in->initialize_from_start();
        while (step != steps.end()) {
            V* states_out = (*step)->forward(*states_in, &num_edmans);
            delete states_in;
            states_in = states_out;
            step++;
        }
        double result = states_in->sum();
        delete states_in;
        return result;
    }

    // This will fit the data the HMM was provided with. It also computes the
    // probability as a side effect, so it returns this in case that is useful
    // to the caller.
    double improve_fit(SequencingModelFitter* fitter) const {
        // There is one less Edman than the number of timesteps, because no
        // Edman is done before the zeroth timestep.
        unsigned int num_edmans = num_timesteps - 1;
        auto step = steps.end();  // const_iterator type
        std::vector<V*> backward_sv;
        backward_sv.reserve(steps.size());
        // For efficiency, backwards_states is in the reverse order of what we
        // would like. Yes this is confusing...
        backward_sv.push_back(create_states_backward());
        backward_sv.back()->initialize_from_finish();
        while (step != steps.begin()) {
            step--;
            V* right_states = backward_sv.back();
            backward_sv.push_back(
                    (*step)->backward(*right_states, &num_edmans));
        }
        double probability = backward_sv.back()->source();
        // We will end up adding NaN results to the fitter if the probability is
        // 0, because the numerators and denominators of parameter estimates
        // will both be 0 - parameter estimates are less than or equal to the
        // full probability always. Then we end up adding 0/0 to the numerator
        // AND the denominator, putting a NaN in both cases.
        //
        // To avoid this, we just return early.
        if (probability == 0.0) {
            return probability;
        }
        auto backward_states = backward_sv.end();  // iterator type
        V* forward_states = create_states_forward();
        forward_states->initialize_from_start();
        while (step != steps.end()) {
            backward_states--;
            (*step)->improve_fit(*forward_states,
                                 **backward_states,
                                 **(backward_states - 1),
                                 num_edmans,
                                 probability,
                                 fitter);
            delete *backward_states;
            V* next_forward_states =
                    (*step)->forward(*forward_states, &num_edmans);
            delete forward_states;
            forward_states = next_forward_states;
            step++;
        }
        delete forward_states;
        delete *(backward_states - 1);
        return probability;
    }
    std::vector<S*> steps;
    unsigned int num_timesteps;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_HMM_GENERIC_HMM_H
