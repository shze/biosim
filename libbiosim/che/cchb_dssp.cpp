#include "che/cchb_dssp.h"
#include <algorithm>

namespace biosim {
  namespace che {
    // ctor from id; id is assumed unique; throws if not found
    cchb_dssp::cchb_dssp(char __id) : _impl(find(__id)) {}
    // ctor from specificity; assumed not unique, constructs first
    cchb_dssp::cchb_dssp(specificity_type __specificity) : _impl(find(__specificity)) {}
    // ctor from id and weight set
    cchb_dssp::cchb_dssp(char __id, weight_map __weights) {
      if(__weights.empty()) { // need at least one element in weight_map
        throw cchb_dssp_data_not_found("weight_map empty");
      } // if

      for(auto const &p : __weights) {
        find(p.first); // just to make sure all are found b/c it throws otherwise
      } // for

      _impl = std::make_shared<data>(__id, __weights);
    } // ctor

    // get identifier
    char cchb_dssp::get_identifier() const { return _impl->_id; }
    // get identifier char
    char cchb_dssp::get_identifier_char() const { return get_identifier(); }
    // get specificity
    cchb_dssp::specificity_type cchb_dssp::get_specificity() const { return _impl->_specificity; }
    // return if an unknown (unspecified) compound; only for unknown, not for gap.
    bool cchb_dssp::is_unknown() const { return _impl->_specificity == specificity_type::unknown; }
    // return if gap
    bool cchb_dssp::is_gap() const { return _impl->_specificity == specificity_type::gap; }
    // get profile weights
    cchb_dssp::weight_map cchb_dssp::get_weights() const { return _impl->_weights; }

    // define ordering
    bool cchb_dssp::operator<(cchb_dssp const &__rhs) const {
      return get_identifier() < __rhs.get_identifier() ||
             (get_identifier() == __rhs.get_identifier() && get_specificity() < __rhs.get_specificity()) ||
             (get_identifier() == __rhs.get_identifier() && get_specificity() == __rhs.get_specificity() &&
              get_weights() < __rhs.get_weights());
    } // operator<
    // equality
    bool cchb_dssp::operator==(cchb_dssp const &__rhs) const { return !(*this < __rhs) && !(__rhs < *this); }

    // (static) get list of all identifiers; static public interface
    std::list<char> cchb_dssp::get_id_list() {
      std::list<char> ids;
      for(auto const &p : data_enum::get_instances()) {
        ids.push_back(p.first[0]); // since we only allow chars we can safely use the first char
      } // for
      return ids;
    } // get_id_list()
    // (static) create a string of identifier chars; static public interface
    std::string cchb_dssp::get_identifier_char_string() {
      std::list<char> id_list(get_id_list());
      return std::string(id_list.begin(), id_list.end());
    } // get_identifier_char_string

    // ctor for primary symbol
    cchb_dssp::data::data(char __id, specificity_type __specificity)
        : _id(__id), _specificity(__specificity), _weights() {
      _weights = {{_id, 1.0}};
    } // ctor
    // ctor for profile symbol
    cchb_dssp::data::data(char __id, weight_map __weights)
        : _id(__id), _specificity(specificity_type::profile), _weights(__weights) {}

    // (static) find compound with id; id is assumed unique; throws if no object was found
    std::shared_ptr<cchb_dssp::data const> cchb_dssp::find(char const &__id) {
      try {
        return data_enum::get_instances().at(std::string(1, __id)).get_object();
      } catch(std::out_of_range &e) {
        throw cchb_dssp_data_not_found("map out of range, no data found with id=" + std::string(1, __id));
      }
    } // find()
    // (static) find compound with specificity; NOT assumed unique, only the first compound is returned
    std::shared_ptr<cchb_dssp::data const> cchb_dssp::find(specificity_type const &__specificity) {
      auto itr = std::find_if(data_enum::get_instances().begin(), data_enum::get_instances().end(),
                              [&](std::pair<std::string, tools::enumerate<std::shared_ptr<data const>>> p) {
        return p.second.get_object()->_specificity == __specificity;
      });

      if(itr == data_enum::get_instances().end()) {
        throw cchb_dssp_data_not_found("no data found with specificity_type=" + std::to_string((int)__specificity));
      } // if

      return itr->second.get_object();
    }
    // (static) constructs data objects of the enumerated instance set; subset of states from DSSP is used, see:
    // http://en.wikipedia.org/wiki/DSSP_%28protein%29
    bool cchb_dssp::initialize() {
      data_enum::add(std::make_shared<data>('H'));
      data_enum::add(std::make_shared<data>('E'));
      data_enum::add(std::make_shared<data>('C', specificity_type::unknown));
      data_enum::add(std::make_shared<data>('-', specificity_type::gap));

      data_enum::add(std::make_shared<data>('G', weight_map({{'H', 1.0}}))); // DSSP 3-10 helix
      data_enum::add(std::make_shared<data>('I', weight_map({{'H', 1.0}}))); // DSSP Pi helix
      data_enum::add(std::make_shared<data>('B', weight_map({{'C', 1.0}}))); // DSSP isolated beta bridge
      data_enum::add(std::make_shared<data>('T', weight_map({{'C', 1.0}}))); // DSSP hydrogen-bonded turn
      data_enum::add(std::make_shared<data>('S', weight_map({{'C', 1.0}}))); // DSSP bend

      return true;
    } // initialize()

    bool cchb_dssp::_initialized(cchb_dssp::initialize()); // initialize static variable
  } // namespace che
} // namespace biosim