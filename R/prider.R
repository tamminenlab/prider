#' Prider
#'
#' @docType package
#' @author Manu Tamminen <mavatam@utu.fi>, Niina Smolander <nijasm@utu.fi>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib prider
#' @name prider
NULL


#' Prepare a primer table for downstream analyses
#'
#' @param input_fasta A string or a DataFrame containing Id and Seq columns.
#' @param primer_length A number specifying length for the designed primers.
#' @param GCcheck A logical; check the GC contents of the primer halves.
#' @param GCsimilarity A number; if GCcheck is performed, this parameter
#'                     determines the maximum proportional GC content
#'                     difference between the primer halves.
#'
#' @return A list containing sequence id conversions, primer matrix
#'         and a list of primers with their target sequences.
#' @importFrom dplyr select
#' @importFrom dplyr as_tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr row_number
#' @importFrom magrittr %>%
#' @importFrom blaster read_fasta
#' @importFrom stringr str_count
prepare_primer_df <- function(input_fasta,
                              primer_length = 20,
                              GCcheck = FALSE,
                              GCsimilarity = 0.1) {
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
    dplyr::select(Seq_no, Id, Seq) %>%
    dplyr::rename(Original_id = Id, Id = Seq_no, Sequence = Seq) %>%
    dplyr::as_tibble(.)

  primers <-
    fasta_table %>%
    dplyr::select(Seq_no, Seq) %>%
    dplyr::rename(Id = Seq_no) %>%
    dplyr::as_tibble(.) %>%
    chunker(window_size = primer_length +1) %>%
    dplyr::as_tibble(.) %>%
    dplyr::select(Seq, Id) %>%
    dplyr::filter(stringr::str_count(Seq, "A|G|C|T") == nchar(Seq))

  if(isTRUE(GCcheck)){
    firsthalf <- substr(primers$Seq, 1, (primer_length/2))
    secondhalf <- substr(primers$Seq, ((primer_length/2)+1), primer_length)
    firsthalf <- str_count(firsthalf, "G|C") / nchar(firsthalf)
    secondhalf <- str_count(secondhalf, "G|C") / nchar(secondhalf)
    dplyr::filter(primers, abs(firsthalf-secondhalf) <= GCsimilarity) -> primers
  }

  primer_df <-
    primers %>%
    dplyr::select(Id, Seq)%>%
    dplyr::group_by(Seq) %>%
    dplyr::summarise(Ids = paste0(sort(Id), collapse=","))

  return(list(conversion_table, primer_df))
}


#' @title new_prider
#' @param x A list
#' @return A prider object
new_prider <- function(x = list()) {
  stopifnot(is.list(x))
  structure(x, class = "prider")
}


#' Prepare a nearly optimal primer coverage for an input FASTA file
#'
#' @title prider
#'
#' @param fasta_file A string. Name of the input FASTA file.
#' @param primer_length A number. Sets the primer length.
#' @param minimum_primer_group_size A number. Sets the minimum number of primers per primer cluster; smaller
#'                                            primer clusters will be discarded.
#' @param minimum_seq_group_size A number. Sets the minimum number of target sequences each primer cluster has to cover.
#' @param cum_cov_decimals A number. Sets the number of decimals for cumulative coverage of primer clusters.
#'                                   More decimals corresponds to more clusters.
#' @param GCcheck A logical; check the GC contents of the primer halves.
#' @param GCsimilarity A number; if GCcheck is performed, this parameter
#'                     determines the maximum proportional GC content
#'                     difference between the primer halves.
#'
#' @return A list containing a sequence conversion table and
#'         a primer coverage table.
#' @examples
#'
#' test_fasta <- system.file("extdata", "test.fasta", package = "prider")
#'
#' primer_designs <- prider(test_fasta)
#'
#' primers(primer_designs)
#'
#' primers(primer_designs)[1]
#'
#' sequences(primer_designs)
#'
#' sequences(primer_designs)[1]
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
#' @importFrom dplyr desc
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
prider <- function(fasta_file,
                   primer_length = 20,
                   minimum_primer_group_size = 10,
                   minimum_seq_group_size = 5,
                   cum_cov_decimals = 2,
                   GCcheck = FALSE,
                   GCsimilarity = 0.1) {

  cat("Preparing primer candidates...\n")
  ag_data <-
    prepare_primer_df(fasta_file, primer_length, GCcheck, GCsimilarity)

  if (!is.character(fasta_file))
    description <- paste0("Primer candidates for DataFrame ", deparse(substitute(fasta_file)), ":\n")
  else
    description <- paste0("Primer candidates for file ", fasta_file, ":\n")

  cat("Clustering primers...\n")
  primer_df_summ <-
    ag_data[[2]] %>%
    dplyr::group_by(Ids) %>%
    dplyr::summarise(Primers = paste0(sort(Seq), collapse=",")) %>%
    dplyr::ungroup()

  abund_clusters <-
    primer_df_summ %>%
    dplyr::mutate(Primer_group_size = lengths(strsplit(primer_df_summ$Primers, ","))) %>%
    dplyr::filter(Primer_group_size >= minimum_primer_group_size)

  if(isTRUE(max(lengths(strsplit(abund_clusters$Ids, ","))) < minimum_seq_group_size)) {
    stop("All sequence group sizes smaller than the minimum_seq_group_size resulting to an empty dataframe.\nPlease make the minimum_seq_group_size parameter smaller")
  }

  cat("Sampling primers...\n")
  primer_draws <-
    abund_clusters %>%
    dplyr::mutate(Primer_group = rownames(abund_clusters)) %>%
    dplyr::mutate(Seq_group_size = lengths(strsplit(abund_clusters$Ids, ","))) %>%
    dplyr::filter(Seq_group_size >= minimum_seq_group_size) %>%
    dplyr::mutate(Sequences = Ids, .keep = "unused")
  primer_draws <- dplyr::arrange(primer_draws, desc(Seq_group_size))

  cat("Eliminating redundancies...\n")
  all_seqs <- sort(unique(unlist(strsplit(primer_draws$Sequences, ","))))
  cum_coverage <- c()
  cum_seqs <- c()
  for (seqs in primer_draws$Sequences) {
    seqs <- unlist(str_split(seqs, ","))
    cum_seqs <- unique(c(cum_seqs, seqs))
    cover <- sum(all_seqs %in% cum_seqs) / length(all_seqs)
    cum_coverage <- c(cum_coverage, cover)
  }

  primer_draws$Cumulative_coverage <- round(cum_coverage, cum_cov_decimals)
  primer_draws <-
    primer_draws %>%
    dplyr::group_by(Cumulative_coverage) %>%
    dplyr::arrange(Cumulative_coverage, desc(Seq_group_size), desc(Primer_group_size), .by_group = TRUE) %>%
    dplyr::slice_max(1) %>%
    dplyr::select(Primer_group, Primer_group_size,
                  Seq_group_size, Cumulative_coverage,
                  Primers, Sequences) %>%
    dplyr::ungroup(.)

  abund_df <-
    primer_draws %>%
    select(Sequences, Primer_group) %>%
    mutate(Sequences = strsplit(Sequences, ",")) %>%
    unnest(Sequences)

  abund_matrix <- table(abund_df$Primer_group, abund_df$Sequences)

  out_matrix <-
    abund_matrix != 0

  if(nrow(out_matrix) > 1 && ncol(out_matrix) > 1){
    out_matrix <-
      out_matrix[,colSums(out_matrix)>0]
    out_matrix <-
      out_matrix[rowSums(out_matrix)>0,]
  } else if (nrow(out_matrix) > 1 && ncol(out_matrix) <= 1){
    out_matrix <-
      out_matrix[rowSums(out_matrix)>0,]
  } else if (nrow(out_matrix) <= 1 && ncol(out_matrix) > 1){
    out_matrix <-
      out_matrix[,colSums(out_matrix)>0]
  }

  excluded_seqs <-
    filter(ag_data[[1]], !(Id %in% colnames(out_matrix)))

  Primer_candidates <-
    filter(primer_draws, primer_draws$Primer_group %in% rownames(out_matrix))

  cat("Done!\n")

  new_prider(
    list(Description = description,
         Conversion_table = ag_data[[1]],
         Primer_candidates = primer_draws,
         Excluded_sequences = excluded_seqs,
         Primer_matrix = out_matrix ))
}

#' @title print.prider
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
  cat("To access the primer candidates,\n")
  cat("please use function primers: primers(.)\n")
  cat("\n")
  cat("To access the sequences covered by the primer candidates,\n")
  cat("please use function sequences: sequences(.)\n")
}


#' @title plot.prider
#' @rdname prider
#' @export
#' @importFrom gplots heatmap.2
plot.prider <- function(prider_obj) {
  matr <- prider_obj$Primer_matrix * 1
  if(ncol(matr) >= 2 && nrow(matr) >= 2){
    gplots::heatmap.2(matr,
                      scale="none",
                      trace="none",
                      col=c("white", "black"),
                      xlab="Sequence Id",
                      ylab="Primer cluster")
  } else {
    cat("Primer_matrix too small to be plotted.")
    cat("Please view the Primer_matrix to see the Ids and clusters.")
  }
}

#' Primers object constructor
#'
#' @title new_primers
#' @param x A tibble
#' @return A primers object
#' @importFrom tibble is_tibble
new_primers <- function(x) {
  stopifnot(tibble::is_tibble(x))
  structure(x, class = "primers")
}


#' Definitions for the S3 methods for the primers classes
#'
#' @title primers
#' @examples
#'
#' test_fasta <- system.file("extdata", "test.fasta", package = "prider")
#'
#' primer_designs <- prider(test_fasta)
#'
#' primers(primer_designs)
#'
#' primers(primer_designs)[1]
#'
#' @export
primers <- function(x) UseMethod("primers")


#' @rdname primers
#' @method primers default
#' @export
primers.default <- function(x, ...)
  warning(paste("Function 'primer' does not know how to handle object of class",
                class(x),
                "and can only be used on class 'prider'."))


#' @rdname primers
#' @method primers prider
#' @export
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
primers.prider <- function(prider_obj) {
  primer_obj <-
    prider_obj$Primer_candidates %>%
    dplyr::select(Primer_group, Seq_group_size, Primer_group_size, Primers) %>%
    dplyr::group_by(Primer_group) %>%
    dplyr::mutate(Primers = str_split(Primers, ",")) %>%
    dplyr::ungroup(.)
  new_primers(primer_obj)
}


#' @rdname primers
#' @export
#' @importFrom dplyr select
print.primers <- function(primer_obj) {
  cat("To access the primers within a group,\n")
  cat("please use subsetting, eg. primers(.)[42]\n")
  cat("\n")
  class(primer_obj) <- "data.frame"
  primer_obj %>%
    dplyr::select(Primer_group, Seq_group_size, Primer_group_size) %>%
    data.frame %>%
    print(.)
}


#' @rdname primers
#' @export
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
`[.primers` <- function(primer_obj, ix) {
  if (is.numeric(ix))
    ix <- as.character(ix)
  class(primer_obj) <- "data.frame"
  dplyr::filter(primer_obj, Primer_group == ix) %>%
    dplyr::pull(Primers) %>%
    .[[1]] %>%
    unique
}


#' Sequences object constructor
#'
#' @title new_sequences
#' @param x A tibble
#' @return A sequences object
#' @importFrom tibble is_tibble
new_sequences <- function(x) {
  stopifnot(tibble::is_tibble(x))
  structure(x, class = "sequences")
}


#' Definitions for the S3 methods for the sequences classes
#'
#' @title sequences
#' @examples
#'
#' test_fasta <- system.file("extdata", "test.fasta", package = "prider")
#'
#' primer_designs <- prider(test_fasta)
#'
#' sequences(primer_designs)
#'
#' sequences(primer_designs)[1]
#'
#' @export
sequences <- function(x) UseMethod("sequences")


#' @rdname sequences
#' @export
sequences.default <- function(x, ...)
  warning(paste("Function 'sequences' does not know how to handle object of class",
                class(x),
                "and can only be used on class 'prider'."))


#' @rdname sequences
#' @export
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
#' @importFrom dplyr left_join
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom stringr str_split
sequences.prider <- function(prider_obj) {
  conversion_tbl <-
    prider_obj$Conversion_table
  sequence_obj <-
    prider_obj$Primer_candidates %>%
    dplyr::select(Primer_group, Seq_group_size, Primer_group_size, Sequences) %>%
    dplyr::group_by(Primer_group) %>%
    dplyr::mutate(Sequences = stringr::str_split(Sequences, ",")) %>%
    dplyr::ungroup(.) %>%
    tidyr::unnest(Sequences) %>%
    dplyr::left_join(conversion_tbl, by=c("Sequences" = "Id")) %>%
    dplyr::select(-Sequences) %>%
    tidyr::nest(Sequences = c(Original_id, Sequence))
  new_sequences(sequence_obj)
}


#' @rdname sequences
#' @export
#' @importFrom dplyr select
print.sequences <- function(sequence_obj) {
  cat("To access the sequences within a group,\n")
  cat("please use subsetting, eg. sequences(.)[42]\n")
  cat("\n")
  class(sequence_obj) <- "data.frame"
  sequence_obj %>%
    dplyr::select(Primer_group, Seq_group_size, Primer_group_size) %>%
    data.frame %>%
    print(.)
}


#' @rdname sequences
#' @export
#' @importFrom dplyr filter
`[.sequences` <- function(sequence_obj, ix) {
  if (is.numeric(ix))
    ix <- as.character(ix)
  class(sequence_obj) <- "data.frame"
  dplyr::filter(sequence_obj, Primer_group == ix) %>%
    .$Sequences %>%
    .[[1]]
}

