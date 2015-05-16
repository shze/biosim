#ifndef math_tensor_h
#define math_tensor_h

#include <vector>
#include <stdexcept>

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
      std::vector<size_t> _vector_increments; // number of elements in each dimension; simplifies finding in _data
      std::vector<T> _data; // stores the actual data

    public:
      // ctor with size of each dimension
      explicit tensor(std::vector<size_t> __sizes) : _sizes(__sizes), _vector_increments(_sizes.size()), _data() {
        size_t total_size(1);
        for(int dim(_sizes.size() - 1); dim >= 0; --dim) {
          total_size *= _sizes[dim];
          _vector_increments[dim] = dim == _sizes.size() - 1 ? 1 : _vector_increments[dim + 1] * _sizes[dim + 1];
        }
        _data = std::vector<T>(total_size);
      }; // ctor

      // returns the rank
      size_t get_rank() const { return _sizes.size(); }
      // returns the size of the specific dimension; can throw the exceptions thrown by _sizes.at()
      size_t get_size(size_t __dimension) { return _sizes.at(__dimension); }
      // returns a const reference to the T at the specified position
      T const &operator()(std::vector<size_t> const &__pos) const {
        if(__pos.size() != _sizes.size()) {
          throw std::out_of_range("tensor_rank=" + std::to_string(get_rank()) + " differs from access_vector_rank=" +
                                  std::to_string(__pos.size()));
        } // if

        size_t real_pos(0);
        for(size_t dim(0); dim < __pos.size(); ++dim) {
          real_pos += __pos[dim] * _vector_increments[dim];
        } // for
        return _data[real_pos];
      } // operator()
      // returns a changeable reference to the T at the specified position
      T &operator()(std::vector<size_t> const &__pos) {
        // avoid code duplication, see http://stackoverflow.com/questions/123758/;  
        return const_cast<T &>(static_cast<const tensor<T> &>(*this).operator()(__pos)); 
      } // operator()
    }; // class tensor
  } // namespace math
} // namespace biosim

#endif // math_tensor_h
