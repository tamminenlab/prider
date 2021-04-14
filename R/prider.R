#' Prider
#' 
#' Prider implements a BLAST-based primer design algorithm for complex sequence sets.
#' 
#' @docType package
#' @author Manu Tamminen <mavatam.@utu.fi>
#' @import Rcpp 
#' @importFrom Rcpp evalCpp
#' @useDynLib prider
#' @name prider
NULL


#' Prepare a primer table for downstream analyses
#'
#' @param input_fasta A string or a DataFrame containing Id and Seq columns
#' @param primer_length A number
#' @return A list containing sequence id conversions, primer matrix
#'         and a list of primers with their target sequences
#' @importFrom dplyr select
#' @importFrom dplyr as_tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr row_number
#' @importFrom magrittr %>%
#' @importFrom blaster read_fasta
prepare_primer_table <- function(input_fasta,
                                 primer_length = 20) {
  if (is.character(input_fasta))
      input_fasta <-
        blaster::read_fasta(input_fasta)

  if (!(all(names(input_fasta) == c("Id", "Seq"))))
    stop("The input must contain Id and Seq columns")

  fasta_table <-
    input_fasta %>%
    dplyr::mutate(Seq_no = paste0("P", dplyr::row_number()))

  conversion_table <-
    fasta_table %>%
    select(Seq_no, Id, Seq) %>%
    dplyr::rename(Original_id = Id, Id = Seq_no, Sequence = Seq) %>%
    dplyr::as_tibble(.)

  primers <-
    fasta_table %>%
    dplyr::select(Seq_no, Seq) %>%
    dplyr::rename(Id = Seq_no) %>%
    dplyr::as_tibble(.) %>%
    chunker(window_size = primer_length) %>%
    dplyr::as_tibble(.) %>%
    dplyr::select(Seq, Id)

  primer_table <-
    primers %>%
    table(.)

  boolean_table <-
    primer_table == 1

  pri_targets <-
    primers %>%
    dplyr::group_by(Seq) %>%
    dplyr::summarise(Seq_group_size = n(),
                     Sequences = paste0(sort(Id), collapse =","))

  return(list(conversion_table, boolean_table, pri_targets))
}


#' Randomly samples primer matrix to recoved the specified coverage
#'
#' @param primer_table A matrix
#' @param coverage A number
#' @return A character vector containing the primer group number ids
sample_coverage <- function(primer_table,
                            coverage = 0.9) {
  row_len <- ncol(primer_table)
  tbl_len <- rownames(primer_table)
  draws <- sample(tbl_len)
  row_acc <- primer_table[draws[1], ]
  acc <- 0
  for (draw in draws[2:length(draws)]) {
    acc <- acc + 1
    row_acc <- row_acc | primer_table[draw, ]
    if (sum(row_acc) / row_len > coverage)
      break
  }
  return(draws[1:acc])
}


#' Prepare a (nearly) optimal primer coverage of the sample set
#'
#' @title prider
#' @param fasta_file A string
#' @param primer_length A number
#' @param coverage A number
#' @param minimum_primer_group_size A number
#' @param minimum_sequence_group_size A number
#' @param draws A number
#' @return A list containing a sequence conversion table and 
#'         a primer coverage table
#' @examples
#' 
#' test_fasta <- system.file("extdata", "test.fasta", package = "prider")
#' 
#' primer_designs <- prider(test_fasta)
#' 
#' @export
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr as_tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr summarise
#' @importFrom dplyr row_number
#' @importFrom dplyr left_join
#' @importFrom dplyr ungroup
#' @importFrom dplyr slice_max
#' @importFrom dplyr arrange
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
prider <- function(fasta_file,
                   primer_length = 20,
                   coverage = 0.9,
                   minimum_primer_group_size = 10,
                   minimum_sequence_group_size = 10,
                   draws = 100) {

  message("Tabulating primers...")
  ag_data <-
    prepare_primer_table(fasta_file, primer_length)

  if (!is.character(fasta_file))
    description <- paste0("Primer candidates for DataFrame ", deparse(substitute(fasta_file)), ":\n")
  else
    description <- paste0("Primer candidates for file ", fasta_file, ":\n")

  message("Clustering primers...")
  big_clusters <-
    ag_data[[3]] %>%
    filter(Seq_group_size >= minimum_sequence_group_size)

  abund_primers <-
    ag_data[[2]][big_clusters$Seq, ]

  primer_clusters <-
    group_primers(abund_primers)

  abund_clusters <-
    primer_clusters %>%
    .$Primer_groups %>%
    dplyr::group_by(Primer_group) %>% 
    dplyr::mutate(Primer_group_size = dplyr::n()) %>% 
    dplyr::filter(Primer_group_size >= minimum_primer_group_size) %>% 
    dplyr::mutate(Primer_group = as.character(Primer_group))

  abund_matrix <-
    primer_clusters %>% 
    .$Primer_matrix %>% 
    .[unique(abund_clusters$Primer_group), ]

  abund_matrix <-
    abund_matrix[, colSums(abund_matrix) != 0]

  message("Sampling primers...")
  primer_draws <-
    purrr::map(1:draws, ~sample_coverage(abund_matrix, coverage)) %>% 
    purrr::map_dfr(~tibble::tibble(Primer_group = .), .id="Draw") %>%
    dplyr::group_by(Draw) %>% 
    dplyr::mutate(Draw_size = dplyr::n()) %>% 
    dplyr::left_join(abund_clusters, by = "Primer_group") %>% 
    dplyr::ungroup(.) %>% 
    dplyr::left_join(ag_data[[3]], by = c("Primer" = "Seq")) %>% 
    dplyr::filter(Draw_size == min(Draw_size)) %>% 
    dplyr::select(-Draw, -Draw_size) %>%
    dplyr::group_by(Primer_group, Primer_group_size, Seq_group_size, Sequences) %>%
    dplyr::summarise(Primers = paste0(sort(Primer), collapse=",")) %>%
    dplyr::ungroup(.) %>% 
    dplyr::arrange(desc(Seq_group_size))

  covered_seqs <-
    primer_draws$Sequences %>%
    unique %>%
    stringr::str_split(",") %>%
    unlist  %>% 
    unique

  excluded_seqs <- 
    filter(ag_data[[1]], !(Id %in% covered_seqs))

  all_seqs <- colnames(abund_matrix)
  cum_coverage <- c()
  cum_seqs <- c()
  for (seqs in primer_draws$Sequences) {
    seqs <- unlist(str_split(seqs, ","))
    cum_seqs <- unique(c(cum_seqs, seqs))
    cover <- sum(all_seqs %in% cum_seqs) / length(all_seqs)
    cum_coverage <- c(cum_coverage, cover)
  }

  primer_draws$Cumulative_coverage <- round(cum_coverage, 2)
  primer_draws <-
    primer_draws %>%
    dplyr::group_by(Cumulative_coverage) %>%
    dplyr::slice_max(1) %>% 
    dplyr::select(Primer_group, Primer_group_size,
           Seq_group_size, Cumulative_coverage,
           Primers, Sequences)

  rows <- unique(primer_draws$Primer_group)
  out_matrix <- abund_matrix[rows, ]

  message("Done!")

  prider_output <- 
    list(Description = description,
         Conversion_table = ag_data[[1]],
         Primer_candidates = primer_draws,
         Excluded_sequences = excluded_seqs,
         Primer_matrix = out_matrix )
  class(prider_output) <- "prider"
  return(prider_output)
}


#' @rdname prider
#' @export
#' @importFrom dplyr count
print.prider <- function(prider_obj) {
  descr <- prider_obj$Description
  input_seqs <- nrow(prider_obj$Conversion)
  excl_seqs <- nrow(prider_obj$Excluded_sequences)
  incl_seqs <- input_seqs - excl_seqs
  primer_candidates <- sum(prider_obj$Primer_candidates$Primer_group_size)
  groups <- nrow(dplyr::count(prider_obj$Primer_candidates, Primer_group))

  total <- paste(input_seqs, "input sequences;", incl_seqs, "included and", excl_seqs, "excluded.\n")
  cands <- paste(primer_candidates, "primer candidates in", groups, "groups.\n")
  cat(descr)
  cat(total)
  cat(cands)
  cat("\n")
  cat("For details:\n")
  cat(".$Primer_candidates\n")
  cat(".$Conversion_table\n")
  cat(".$Excluded_sequences\n")
}


#' @rdname prider
#' @export
#' @importFrom gplots heatmap.2
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 coord_flip
#' @importFrom dplyr select
#' @importFrom tidyr pivot_longer
plot.prider <- function(prider_obj, plot_type = "heatmap") {
  if (plot_type == "heatmap") {
    matr <- prider_obj$Primer_matrix * 1
    gplots::heatmap.2(matr,
                      scale="none",
                      trace="none",
                      col=c("white", "black"),
                      xlab="Sequence Id",
                      ylab="Primer cluster")
  } else if (plot_type == "barchart") {
    prider_obj$Primer_candidates %>%
      dplyr::select(Primer_group, Primer_group_size, Seq_group_size) %>%
      tidyr::pivot_longer(cols=-"Primer_group", names_to="Var", values_to="Val") %>%
      ggplot2::ggplot(ggplot2::aes(x=Primer_group, y=Val, fill=Var)) +
      ggplot2::geom_bar(stat="identity", position="dodge") +
      ggplot2::coord_flip()
  } else {
    stop("Currently only heatmap and barchart are supported!")
  }
}

