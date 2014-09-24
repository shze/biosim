#include "che/cchb.h"
#include <algorithm>

namespace biosim {
  namespace che {
    // construct cc with id; id is assumed unique; throws if not found
    cchb::cchb(std::string __id) : _impl(find(__id)) {}
    // ctor from specificity ; assumed not unique, constructs first
    cchb::cchb(specificity_type __specificity) : _impl(find(__specificity)) {}
    // ctor from id and weight set
    cchb::cchb(std::string __id, weight_map __weights) {
      if(__weights.empty()) { // need at least one element in weight_map
        throw cchb_data_not_found("weight_map empty");
      } // if

      for(auto const &p : __weights) {
        find(p.first); // just to make sure all are found b/c it throws otherwise
      } // for

      _impl = std::make_shared<data>(__id, __weights);
    } // ctor

    // get identifier
    std::string cchb::get_identifier() const { return _impl->_id; }
    // get specificity
    cchb::specificity_type cchb::get_specificity() const { return _impl->_specificity; }
    // return if an unknown (unspecified) compound; only for unknown, not for gap.
    bool cchb::is_unknown() const { return _impl->_specificity == specificity_unknown; }
    // return if gap
    bool cchb::is_gap() const { return _impl->_specificity == specificity_gap; }
    // get profile weights
    cchb::weight_map cchb::get_weights() const { return _impl->_weights; }

    // define ordering
    bool cchb::operator<(cchb const &__rhs) const {
      return get_identifier() < __rhs.get_identifier() ||
             (get_identifier() == __rhs.get_identifier() && get_specificity() < __rhs.get_specificity()) ||
             (get_identifier() == __rhs.get_identifier() && get_specificity() == __rhs.get_specificity() &&
              get_weights() < __rhs.get_weights());
    } // operator<
    // equality
    bool cchb::operator==(cchb const &__rhs) const { return !(*this < __rhs) && !(__rhs < *this); }

    // (static) get list of all identifiers; static public interface
    std::list<std::string> cchb::get_id_list() {
      std::list<std::string> ids;
      for(auto p : data_enum::get_instances()) {
        ids.push_back(p.first);
      } // for
      return ids;
    } // get_id_list()

    // ctor for primary symbol
    cchb::data::data(std::string __id, specificity_type __specificity)
        : _id(__id), _specificity(__specificity), _weights() {
      _weights = {{_id, 1.0}};
    } // ctor
    // ctor for profile symbol
    cchb::data::data(std::string __id, weight_map __weights)
        : _id(__id), _specificity(specificity_defined), _weights(__weights) {}

    // (static) find compound with id; id is assumed unique; throws if no object was found
    std::shared_ptr<cchb::data const> cchb::find(std::string __id) {
      try {
        return data_enum::get_instances().at(__id).get_object();
      } catch(std::out_of_range &e) {
        throw cchb_data_not_found("map out of range, no data found with id=" + __id);
      }
    } // find()
    // (static) find compound with specificity; NOT assumed unique, only the first compound is returned
    std::shared_ptr<cchb::data const> cchb::find(specificity_type __specificity) {
      auto itr = std::find_if(data_enum::get_instances().begin(), data_enum::get_instances().end(),
                              [&](std::pair<std::string, tools::enumerate<std::shared_ptr<data const>>> p) {
        return p.second.get_object()->_specificity == __specificity;
      });

      if(itr == data_enum::get_instances().end()) {
        throw cchb_data_not_found("no data found with specificity_type=" + std::to_string(__specificity));
      } // if

      return itr->second.get_object();
    }
    // (static) constructs data objects of the enumerated instance set; subset of states from DSSP is used, see:
    // http://en.wikipedia.org/wiki/DSSP_%28protein%29
    bool cchb::initialize() {
      data_enum::add(std::make_shared<data>("H"));
      data_enum::add(std::make_shared<data>("E"));
      data_enum::add(std::make_shared<data>("C", specificity_unknown));
      data_enum::add(std::make_shared<data>("-", specificity_gap));

      return true;
    } // initialize()

    bool cchb::_initialized(cchb::initialize()); // initialize static variable
  } // namespace che
} // namespace biosim