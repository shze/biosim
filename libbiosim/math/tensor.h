#ifndef math_tensor_h
#define math_tensor_h

#include <vector>
#include <stdexcept>
#include "tools/incrementor.h"

namespace biosim {
  namespace math {
    // tensor class for arbitary rank, but rank and size is fixed on construction;
    // for the math, see http://math.stackexchange.com/questions/73144/;
    // design:
    // set rank and sizes of dimensions as parameters of ctor; do not use variadic function args, b/c their number
    // cannot be determined, risking overrunning the end of the list; do not use a rank template parameter, b/c this
    // avoids assigning tensors of different ranks which is not a problem; do not use boost::multi_array b/c it allows
    // access to the subdimensions, which is not desired here.
    template <class T>
    class tensor {
    private:
      std::vector<size_t> _sizes; // size of each dimension
      std::vector<size_t> _strides; // number of elements in each dimension; simplifies finding in _data
      std::vector<T> _data; // stores the actual data

    public:
      // ctor with size of each dimension
      explicit tensor(std::vector<size_t> __sizes) : _sizes(__sizes), _strides(_sizes.size()), _data() {
        size_t total_size(1);
        for(int dim(_sizes.size() - 1); dim >= 0; --dim) {
          total_size *= _sizes[dim];
          _strides[dim] = dim == _sizes.size() - 1 ? 1 : _strides[dim + 1] * _sizes[dim + 1];
        }
        _data = std::vector<T>(total_size);
      }; // ctor

      // returns the rank
      size_t get_rank() const { return _sizes.size(); }
      // returns the size of the specific dimension; can throw the exceptions thrown by _sizes.at()
      size_t get_size(size_t __dimension) { return _sizes.at(__dimension); }
      // returns a const reference to the T at the specified position
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
      // returns the subtensor with the given dimensions (0..rank - 1) at the given positions; elements are copied
      tensor<T> sub(std::vector<size_t> const &__subt_dimensions, std::vector<size_t> const &__subt_positions) {
        DEBUG << "Creating subtensor with rank=" << std::to_string(__subt_dimensions.size())
              << " and positions=" << std::to_string(__subt_positions.size())
              << " from tensor with rank=" << get_rank();
        if(__subt_dimensions.size() + __subt_positions.size() != get_rank()) {
          throw std::out_of_range("cannot get subtensor with rank=" + std::to_string(__subt_dimensions.size()) +
                                  " and positions=" + std::to_string(__subt_positions.size()) +
                                  " from tensor with rank=" + std::to_string(get_rank()));
        } // if

        std::vector<size_t> subt_sizes; // create and fill vector with sizes for subtensor
        std::for_each(__subt_dimensions.begin(), __subt_dimensions.end(),
                      [&](size_t const &__dim) { subt_sizes.push_back(get_size(__dim)); });
        tensor<T> subt(subt_sizes); // create subtensor

        if(__subt_dimensions.size() == 0) { // incrementor throws with empty alphabets for rank=0, so shortcut here
          DEBUG << "Directly returning subtensor with rank=0";
          subt({}) = this->operator()(__subt_positions);
          return subt;
        } // if

        // create alphabets and start position for incrementor
        std::vector<std::vector<size_t>> t_alphabets, subt_alphabets;
        std::vector<size_t> t_pos(get_rank(), 0), subt_pos(__subt_dimensions.size(), 0);
        // do not test current_subtensor_position in the condition, b/c the loop needs to continue with dim even if it
        // is at the end; it never runs past the end b/c of the test/throw at the beginning
        for(size_t dim(0), current_subt_pos(0); dim < get_rank(); ++dim) {
          if(std::find(__subt_dimensions.begin(), __subt_dimensions.end(), dim) == __subt_dimensions.end()) {
            DEBUG << "Add fixed alphabet from subtensor_position=" << current_subt_pos
                  << " for tensor_dimension=" << dim;
            t_alphabets.push_back(std::vector<size_t>({__subt_positions[current_subt_pos]}));
            t_pos[dim] = __subt_positions[current_subt_pos];
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
        tools::incrementor<std::vector<size_t>> tensor_inc(t_alphabets), subtensor_inc(subt_alphabets);
        bool done(false);
        while(!done) {
          subt(subt_pos) = this->operator()(t_pos);

          try {
            DEBUG << "Incrementing tensor position";
            t_pos = tensor_inc.next(t_pos);
            DEBUG << "Incrementing subtensor position";
            subt_pos = subtensor_inc.next(subt_pos);
          } // try
          catch(std::overflow_error &e) {
            done = true;
          } // catch
        } // while

        return subt;
      } // sub()
    }; // class tensor
  } // namespace math
} // namespace biosim

#endif // math_tensor_h
