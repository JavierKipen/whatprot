/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_TENSOR_H
#define WHATPROT_TENSOR_TENSOR_H

// Standard C++ library headers:
#include <initializer_list>

// Local project headers:
#include "tensor/const-tensor-iterator.h"
#include "tensor/const-tensor-vector-iterator.h"
#include "tensor/tensor-iterator.h"
#include "tensor/tensor-vector-iterator.h"
#include "util/kd-range.h"

namespace whatprot {

class Tensor {
public:
    Tensor(unsigned int order, const unsigned int* shape);
    Tensor(const KDRange& range);
    Tensor(Tensor&& other);
    ~Tensor();
    double& operator[](const unsigned int* loc);
    // This next function is probably not useful in production, but is very
    // helpful for readable tests. Note that, unfortunately, lines of code like
    //   - BOOST_TEST(t[{1, 2}] == 314);
    // will not work because they mess up the Boost testing templates. Instead
    // do something like
    //   - BOOST_TEST((t[{1, 2}]) == 314);
    // Notice the additional set of parenthesis in the corrected example line.
    double& operator[](std::initializer_list<unsigned int> loc);
    TensorIterator* iterator(const KDRange& range);
    ConstTensorIterator* const_iterator(const KDRange& range) const;
    TensorVectorIterator* vector_iterator(const KDRange& range,
                                          unsigned int vector_dimension);
    ConstTensorVectorIterator* const_vector_iterator(
            const KDRange& range, unsigned int vector_dimension) const;
    double sum() const;
    double sum(const KDRange& range) const;

    double* values;
    KDRange range;
    int* strides;
    unsigned int size;
    unsigned int order;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_TENSOR_H
