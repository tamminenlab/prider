#include <Rcpp.h>
#include <string>
#include <sstream>
#include <set>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
DataFrame unique_network(DataFrame network)
{
  CharacterVector left = network[0];
  CharacterVector right = network[1];

  std::unordered_set< std::string > netw_set;

  for (int i = 0; i < left.size(); i++) {
    std::string left_entry { left[i] };
    std::string right_entry { right[i] };
    std::string netw_entry {left_entry + "," + right_entry};
    netw_set.insert(netw_entry);
  }

  std::vector < std::string > left_nodes;
  std::vector < std::string > right_nodes;

  for (std::string pair : netw_set)
    {
      std::string left_node { pair.substr(0, pair.find(",")) };
      std::string right_node { pair.substr(pair.find(",") + 1, pair.size()) };
      left_nodes.push_back(left_node);
      right_nodes.push_back(right_node);
    }

  DataFrame seq_table = DataFrame::create(Named("Left") = left_nodes,
                                          Named("Right") = right_nodes);

  return seq_table;

}
