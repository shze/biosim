#ifndef math_tensor_h
#define math_tensor_h

#include <vector>
#include <stdexcept>
#include <iterator>
#include "tools/mapper.h"
#include "tools/log.h"
#include "tools/format.h"

namespace biosim {
  namespace math {
    // tensor class for arbitary rank, but rank and size is fixed on construction;
    // for the math, see http://math.stackexchange.com/questions/73144/;
    // design:
    // set rank and sizes of dimensions as parameters of ctor; do not use variadic function args, b/c their number
    // cannot be determined, risking overrunning the end of the list; do not use a rank template parameter, b/c this
    // avoids assigning tensors of different ranks which is not a problem; do not use boost::multi_array b/c it allows
    // access to the subdimensions, which is not desired here;
    // sizes are ordered from lowest to highest dimension, i.e. a tensor(n_0, n_1) of rank 2 (a matrix) looks like:
    //  a(0, 0)       a(1, 0)       ... a(n_0 - 1, 0)
    //  a(0, 1)       a(1, 1)       ... a(n_0 - 1, 1)
    //  ...           ...           ... ...
    //  a(0, n_1 - 1) a(1, n_1 - 1) ... a(n_0 - 1, n_1 - 1)
    // that is a(i_0, i_1) with i_0 = 0..n_0 - 1 and i_1 = 0..n_1 - 1
    // and a tensor(n_0, n_1, n_2) of rank 3 looks like a(i_0, i_1, i_2) with i_k = 0..n_k - 1;
    // this is the same dimension order as boost matrix uses.
    template <class T>
    class tensor {
    public:
      // ctor with size of each dimension
      explicit tensor(std::vector<size_t> __sizes) : _sizes(__sizes), _strides(_sizes.size()), _data() {
        size_t total_size(1);
        for(size_t dim(0); dim < _sizes.size(); ++dim) { // first to last, i.e. lowest to highest dimension
          total_size *= _sizes[dim];
          _strides[dim] = dim == 0 ? 1 : _strides[dim - 1] * _sizes[dim - 1];
        }
        _data = std::vector<T>(total_size);
      }; // ctor

      // returns the rank
      size_t get_rank() const { return _sizes.size(); }
      // returns the size of the given dimension (0 = lowest, rank-1 = highest), 0 if tensor doesn't have this dimension
      size_t get_size(size_t __dimension) const { return __dimension < get_rank() ? _sizes[__dimension] : 0; }
      // converts dimension sizes into alphabets that can be used to iterate using mapper
      tools::mapper<std::vector<size_t>>::alphabet_container get_mapper_alphabets() const {
        tools::mapper<std::vector<size_t>>::alphabet_container alphabets;
        for(auto b : _sizes) {
          std::vector<size_t> alphabet(b);
          size_t element(0);
          std::generate(alphabet.begin(), alphabet.end(), [&] { return element++; });
          alphabets.push_back(alphabet);
        } // for
        return alphabets;
      } // get_mapper_alphabets()
      // returns a const reference to the T at the given position (n_0, .., n_{rank-1})
      T const &operator()(std::vector<size_t> const &__pos) const {
        if(__pos.size() != get_rank()) {
          throw std::out_of_range("tensor_rank=" + std::to_string(get_rank()) + " differs from access_vector_rank=" +
                                  std::to_string(__pos.size()));
        } // if
        // check each dim is within size limits, only checking real_pos is insufficient
        for(size_t dim(0); dim < get_rank(); ++dim) {
          if(__pos[dim] >= _sizes[dim]) {
            throw std::out_of_range("cannot access position>=size with position=" + std::to_string(__pos[dim]) +
                                    " and size=" + std::to_string(_sizes[dim]) + " for dimension=" +
                                    std::to_string(dim));
          } // if
        } // for

        size_t real_pos(0);
        for(size_t dim(0); dim < __pos.size(); ++dim) {
          real_pos += __pos[dim] * _strides[dim];
        } // for
        return _data[real_pos];
      } // operator()
      // returns a changeable reference to the T at the specified position
      T &operator()(std::vector<size_t> const &__pos) {
        // avoid code duplication, see http://stackoverflow.com/questions/123758/;
        return const_cast<T &>(static_cast<const tensor<T> &>(*this).operator()(__pos));
      } // operator()
      // returns the subtensor with the given list of dimensions at the given positions; elements are copied;
      // subt_dimensions does not need to be sorted, it is sorted before use
      tensor<T> sub(std::vector<size_t> __subt_dimensions, std::vector<size_t> const &__subt_positions) const {
        DEBUG << "Creating subtensor with rank=" << std::to_string(__subt_dimensions.size())
              << " and positions=" << std::to_string(__subt_positions.size())
              << " from tensor with rank=" << get_rank();
        if(__subt_dimensions.size() + __subt_positions.size() != get_rank()) {
          throw std::out_of_range("cannot get subtensor with rank=" + std::to_string(__subt_dimensions.size()) +
                                  " and positions=" + std::to_string(__subt_positions.size()) +
                                  " from tensor with rank=" + std::to_string(get_rank()));
        } // if

        std::sort(__subt_dimensions.begin(), __subt_dimensions.end()); // sort from lowest to highest
        std::vector<size_t> subt_sizes; // create and fill vector with sizes for subtensor
        std::for_each(__subt_dimensions.begin(), __subt_dimensions.end(),
                      [&](size_t const &__dim) { subt_sizes.push_back(get_size(__dim)); }); // insert in order
        tensor<T> subt(subt_sizes); // create subtensor

        // create alphabets for mapper
        std::vector<std::vector<size_t>> t_alphabets, subt_alphabets;
        // do not test current_subtensor_position in the condition, b/c the loop needs to continue with dim even if it
        // is at the end; it never runs past the end b/c of the test/throw at the beginning
        for(size_t dim(0), current_subt_pos(0); dim < get_rank(); ++dim) {
          if(std::find(__subt_dimensions.begin(), __subt_dimensions.end(), dim) == __subt_dimensions.end()) {
            DEBUG << "Add fixed alphabet from subtensor_position=" << current_subt_pos
                  << " for tensor_dimension=" << dim;
            t_alphabets.push_back(std::vector<size_t>({__subt_positions[current_subt_pos]}));
            ++current_subt_pos; // only increment if a pos is used, i.e. if this dim is not part of the subtensor
          } // if
          else {
            DEBUG << "Add full alphabet for tensor_dimension=" << dim;
            std::vector<size_t> tensor_alphabet_full(get_size(dim));
            size_t element(0);
            std::generate(tensor_alphabet_full.begin(), tensor_alphabet_full.end(), [&] { return element++; });
            t_alphabets.push_back(tensor_alphabet_full); // add alphabet to tensor
            subt_alphabets.push_back(tensor_alphabet_full); // and subtensor
          } // else
        } // for

        // iterate/increment and copy
        tools::mapper<std::vector<size_t>> tensor_m(t_alphabets), subtensor_m(subt_alphabets);
        std::vector<size_t> t_pos, subt_pos;
        for(size_t t_i(tensor_m.get_min()), subt_i(subtensor_m.get_min());
            t_i <= tensor_m.get_max() && subt_i <= subtensor_m.get_max(); ++t_i, ++subt_i) {
          t_pos = tensor_m.encode(t_i);
          subt_pos = subtensor_m.encode(subt_i);
          DEBUG << "Mapping tensor_position->subtensor_position: (" << tools::to_string(t_pos) << ")->("
                << tools::to_string(subt_pos) << ")";
          subt(subt_pos) = this->operator()(t_pos);
        } // for

        return subt;
      } // sub()

    private:
      std::vector<size_t> _sizes; // size of each dimension
      std::vector<size_t> _strides; // number of elements in each dimension; simplifies finding in _data
      std::vector<T> _data; // stores the actual data
    }; // class tensor

    // output operator for tensor; inline to avoid multiple definitions error when .h is included multiple times
    template <class T>
    inline std::ostream &operator<<(std::ostream &__out, tensor<T> const &__tensor) {
      __out << "tensor: rank=" << __tensor.get_rank() << ", sizes=(";
      if(__tensor.get_rank() > 0) {
        __out << __tensor.get_size(0);
        for(size_t dim(1); dim < __tensor.get_rank(); ++dim) {
          __out << ", " << __tensor.get_size(dim);
        } // for
      } // if
      __out << ")\n";

      // special cases: rank {0, 1}
      if(__tensor.get_rank() == 0) {
        __out << __tensor({});
        return __out;
      } // if
      if(__tensor.get_rank() == 1) {
        for(size_t pos(0); pos < __tensor.get_size(0); ++pos) {
          __out << __tensor({pos}) << " ";
        } // for
        return __out;
      } // if

      // at this point rank is >= 2: prepare for creating and iterating over all rank 2 subtensors
      std::vector<size_t> subt_dim{0, 1}; // create subtensor dimension vector: iterate over lowest two, i.e. 0 and 1
      // iterate over all subtensors with rank 2 by creating a mapper to iterate over all dims >= 2, i.e. begin+2..end
      tools::mapper<std::vector<size_t>>::alphabet_container alphabets(__tensor.get_mapper_alphabets());
      tools::mapper<std::vector<size_t>>::alphabet_container reduced_alphabets(alphabets.begin() + 2, alphabets.end());
      tools::mapper<std::vector<size_t>> m(reduced_alphabets);

      for(size_t i(m.get_min()); i <= m.get_max(); ++i) {
        std::vector<size_t> subt_pos(m.encode(i));
        if(subt_pos.size() > 0) {
          __out << "subtensor(x, y";
          std::for_each(subt_pos.begin(), subt_pos.end(), [&](size_t const &__dim) { __out << ", " << __dim; });
          __out << ")\n";
        } // if
        tensor<T> subt(__tensor.sub(subt_dim, subt_pos));
        for(size_t pos1(0); pos1 < subt.get_size(1); ++pos1) {
          for(size_t pos0(0); pos0 < subt.get_size(0); ++pos0) {
            __out << subt({pos0, pos1}) << " ";
          } // for
          __out << "\n";
        } // for
      } // for

      return __out;
    } // operator<<()
  } // namespace math
} // namespace biosim

#endif // math_tensor_h
