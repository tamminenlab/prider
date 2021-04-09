#include <map>
#include <string>
#include <Rcpp.h>

using namespace Rcpp;

std::string compress_vector(LogicalVector vect) {
  std::string acc;
  for (bool element : vect)
    if (element)
      acc += "1";
    else
      acc += "0";
  return acc;
}

// [[Rcpp::export]]
List group_primers(LogicalMatrix primer_tbl) {

  CharacterVector matr_colnames = colnames(primer_tbl);
  std::unordered_map<std::string, std::vector<std::string>> primer_acc;
  std::unordered_map<std::string, LogicalVector> row_acc;
  std::vector<std::string> rnames = as<std::vector<std::string>>(rownames(primer_tbl));
  int row_len;

  for (int ix = 0; ix < rnames.size(); ix++) {
    LogicalVector row = primer_tbl( ix, _ );
    std::string row_str = compress_vector(row);
    std::string primer = rnames[ix];
    primer_acc[row_str].push_back(primer);
    row_acc[row_str] = row;
    row_len = row.size();
  }

  std::vector<std::string> primers;
  std::vector<int> ixs;
  std::vector<int> uniq_ixs;
  LogicalMatrix logmat(primer_acc.size(), row_len);
  int ix { 0 };

  for (auto elem : primer_acc) {
    std::vector<std::string> primer = elem.second;
    for (std::string pri_elem : primer) {
      ixs.push_back(ix);
      primers.push_back(pri_elem);
    }
    LogicalVector row = row_acc[elem.first];
    uniq_ixs.push_back(ix);
    logmat(ix, _) = row;
    ix++;
  }

  NumericVector rnames_char = wrap(uniq_ixs);
  rownames(logmat) = rnames_char;
  colnames(logmat) = matr_colnames;

  DataFrame primer_groups = DataFrame::create(Named("Primer_group") = ixs,
                                              Named("Primer") = primers);

  return(List::create(Named("Primer_matrix") = logmat,
                      Named("Primer_groups") = primer_groups));
}
