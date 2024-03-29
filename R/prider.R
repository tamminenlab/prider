#' Prider
#'
#' @docType package
#' @author Manu Tamminen <mavatam@utu.fi>, Niina Smolander <nijasm@utu.fi>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib prider
#' @name prider
NULL

utils::globalVariables(c("Primer_group", "Primers", ".", "Seq_no", "Id", "Seq", "Ids",
    "Primer_group_size", "Seq_group_size", "Cumulative_coverage", "Sequences", "Original_id",
    "x", "Sequence"))

#' Prepare a primer table for downstream analyses
#'
#' @param input_fasta A string. Name or filepath of the input FASTA file.
#' @param primer_length A number. Sets the primer length. For applications
#'        involving two adjacent probes, the value should be set to two-fold
#'        the length of a single probe.
#' @param GCcheck A logical. If TRUE, checks the GC contents of the primers and
#'        filters based on GCmin and GCmax.
#' @param GCmin A decimal. If GCcheck is performed, this parameter determines
#'        the minimum proportional GC content.
#' @param GCmax A decimal. If GCcheck is performed, this parameter determines
#'        the maximum proportional GC content.
#' @param GChalves A logical. If TRUE, checks the GC contents separately for
#'        both halves of the primers and filters based on GCsimilarity. For
#'        example for applications involving two adjacent probes.
#' @param GCsimilarity A number. If GChalves is performed, this parameter
#'                     determines the maximum proportional GC content
#'                     difference between the primer halves.
#' @param NTkeep A string. Filters the primers based on the nucleotides. \cr
#'        "basic" = keeps only the primers with G, C, T or A. \cr
#'        "N" = keeps primers with G, C, T, A and N. \cr
#'        "any" = keeps primers with the IUPAC nucleotide code characters. \cr
#'        "all" = keeps all primers.
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
prepare_primer_df <- function(input_fasta, primer_length = 20, NTkeep = "basic", GCcheck = FALSE, GCmin = 0.4,
    GCmax = 0.6, GChalves = FALSE, GCsimilarity = 0.1) {
    if (is.character(input_fasta))
        input_fasta <- blaster::read_fasta(input_fasta, non_standard_chars = "ignore")

    if (!(all(names(input_fasta) == c("Id", "Seq"))))
        stop("The input must contain Id and Seq columns")

    fasta_table <- input_fasta %>%
        dplyr::mutate(Seq_no = paste0("P", dplyr::row_number()))

    conversion_table <- fasta_table %>%
        dplyr::select(Seq_no, Id, Seq) %>%
        dplyr::rename(Original_id = Id, Id = Seq_no, Sequence = Seq) %>%
        dplyr::as_tibble(.)

    primers <- fasta_table %>%
        dplyr::select(Seq_no, Seq) %>%
        dplyr::rename(Id = Seq_no) %>%
        dplyr::as_tibble(.) %>%
        chunker(window_size = primer_length + 1) %>%
        dplyr::as_tibble(.) %>%
        dplyr::select(Seq, Id)
    
    switch(NTkeep,
           "basic" = primers <- primers %>% 
             dplyr::filter(stringr::str_count(Seq, "A|G|C|T") == nchar(Seq)),
           "N" = primers <- primers %>% 
             dplyr::filter(stringr::str_count(Seq, "A|G|C|T|N") == nchar(Seq)),
           "any" = primers <- primers %>% 
             dplyr::filter(stringr::str_count(Seq, "A|G|C|T|N|R|Y|S|W|K|M|B|D|H|V|U") == nchar(Seq)),
           "all" = primers <- primers
    )

    if (isTRUE(GCcheck)) {
        content <- stringr::str_count(primers$Seq, "G|C")/nchar(primers$Seq)
        dplyr::filter(primers, content >= GCmin & content <= GCmax) -> primers
    }

    if (isTRUE(GChalves)) {
        firsthalf <- substr(primers$Seq, 1, (primer_length/2))
        secondhalf <- substr(primers$Seq, ((primer_length/2) + 1), primer_length)
        firsthalf <- stringr::str_count(firsthalf, "G|C")/nchar(firsthalf)
        secondhalf <- stringr::str_count(secondhalf, "G|C")/nchar(secondhalf)
        dplyr::filter(primers, abs(firsthalf - secondhalf) <= GCsimilarity) ->
            primers
    }

    primer_df <- primers %>%
        dplyr::select(Id, Seq) %>%
        dplyr::group_by(Seq) %>%
        dplyr::summarise(Ids = paste0(sort(Id), collapse = ","))

    return(list(conversion_table, primer_df))
}


#' @title new_prider
#' @param x A list
#' @return A prider object
new_prider <- function(x = list()) {
    stopifnot(is.list(x))
    structure(x, class = "prider")
}


#' Prepare a nearly optimal primer coverage for an input FASTA file.
#'
#' @title prider
#'
#' @param fasta_file A string. Name or filepath of the input FASTA file.
#' @param primer_length A number. Sets the primer length. For applications
#'        involving two adjacent probes, the value should be set to two-fold
#'        the length of a single probe.
#' @param minimum_primer_group_size A number. Sets the minimum number of primers
#'        per primer cluster; smaller primer clusters will be discarded.
#' @param minimum_seq_group_size A number. Sets the minimum number of sequences
#'        each primer cluster has to cover.
#' @param cum_cov_decimals A number. Sets the number of decimals for cumulative
#'        coverage of primer clusters. Generally, lower value corresponds to
#'        less clusters and higher value to more clusters in the output. If
#'        the clusters do not cover the input sequences sufficiently, increasing
#'        this value may increase the coverage. If the clusters overlap too
#'        much, lowering the value may reduce this effect. Recommended range 1-4.
#' @param GCcheck A logical. If TRUE, checks the GC contents of the primers and
#'        filters based on GCmin and GCmax.
#' @param GCmin A decimal. If GCcheck is performed, this parameter determines
#'        the minimum proportional GC content.
#' @param GCmax A decimal. If GCcheck is performed, this parameter determines
#'        the maximum proportional GC content.
#' @param GChalves A logical. If TRUE, checks the GC contents separately for
#'        both halves of the primers and filters based on GCsimilarity. Used
#'        for example for applications involving two adjacent probes.
#' @param GCsimilarity A decimal. If GChalves is performed, this parameter
#'                     determines the maximum proportional GC content
#'                     difference between the primer halves.
#' @param NTkeep A string. Filters the primers based on the nucleotides. \cr
#'        "basic" = keeps only the primers with G, C, T or A. \cr
#'        "N" = keeps primers with G, C, T, A and N. \cr
#'        "any" = keeps primers with the IUPAC nucleotide code characters. \cr
#'        "all" = keeps all primers.
#'
#' @return A list containing a sequence conversion table, primer candidates table,
#'         excluded sequences table and a primer coverage table.
#' @examples
#' test_fasta <- system.file('extdata', 'test.fasta', package = 'prider')
#'
#' # Runs Prider with the default values:
#' primer_designs <- prider(test_fasta)
#'
#' # Returns all the primers:
#' primers(primer_designs)
#' # Returns the primers of a specific primer group:
#' primers(primer_designs)[1]
#' # Returns all the sequences:
#' sequences(primer_designs)
#' # Returns the sequence of a specific Id:
#' sequences(primer_designs)[1]
#' # Plots the primers groups and the target sequences as a heatmap:
#' plot(primer_designs)
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
#' @importFrom dplyr slice_head
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
prider <- function(fasta_file, primer_length = 20, minimum_primer_group_size = 10,
    minimum_seq_group_size = 2, cum_cov_decimals = 2, NTkeep = "basic", GCcheck = FALSE, GCmin = 0.4,
    GCmax = 0.6, GChalves = FALSE, GCsimilarity = 0.1) {

    message("Preparing primer candidates...\n")
    ag_data <- prepare_primer_df(fasta_file, primer_length, NTkeep, GCcheck, GCmin, GCmax,
        GChalves, GCsimilarity)

    if (!is.character(fasta_file))
        description <- paste0("Primer candidates for DataFrame ", deparse(substitute(fasta_file)),
            ":\n") else description <- paste0("Primer candidates for file ", fasta_file, ":\n")

    message("Clustering primers...\n")
    primer_df_summ <- ag_data[[2]] %>%
        dplyr::group_by(Ids) %>%
        dplyr::summarise(Primers = paste0(sort(Seq), collapse = ",")) %>%
        dplyr::ungroup()

    abund_clusters <- primer_df_summ %>%
        dplyr::mutate(Primer_group_size = lengths(strsplit(primer_df_summ$Primers,
            ","))) %>%
        dplyr::filter(Primer_group_size >= minimum_primer_group_size)

    if (isTRUE(nrow(abund_clusters) < 1)) {
        stop("No primer candidates left in the data after clustering. \nPlease change the parameters or check the input FASTA file.")
    }

    if (isTRUE(max(lengths(strsplit(abund_clusters$Ids, ","))) < minimum_seq_group_size)) {
        stop("All sequence group sizes smaller than the minimum_seq_group_size resulting to an empty dataframe.\nPlease make the minimum_seq_group_size parameter smaller or change other parameters.")
    }

    message("Sampling primers...\n")
    primer_draws <- abund_clusters %>%
        dplyr::mutate(Primer_group = rownames(abund_clusters)) %>%
        dplyr::mutate(Seq_group_size = lengths(strsplit(abund_clusters$Ids, ","))) %>%
        dplyr::filter(Seq_group_size >= minimum_seq_group_size) %>%
        dplyr::mutate(Sequences = Ids, .keep = "unused")
    primer_draws <- dplyr::arrange(primer_draws, desc(Seq_group_size))

    message("Eliminating redundancies...\n")
    all_seqs <- sort(unique(unlist(strsplit(primer_draws$Sequences, ","))))
    cum_coverage <- c()
    cum_seqs <- c()
    for (seqs in primer_draws$Sequences) {
        seqs <- unlist(stringr::str_split(seqs, ","))
        cum_seqs <- unique(c(cum_seqs, seqs))
        cover <- sum(all_seqs %in% cum_seqs)/length(all_seqs)
        cum_coverage <- c(cum_coverage, cover)
    }

    primer_draws$Cumulative_coverage <- round(cum_coverage, cum_cov_decimals)
    primer_draws <- primer_draws %>%
        dplyr::group_by(Cumulative_coverage) %>%
        dplyr::arrange(Cumulative_coverage, desc(Seq_group_size), desc(Primer_group_size),
            .by_group = TRUE) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::select(Primer_group, Primer_group_size, Seq_group_size, Cumulative_coverage,
            Primers, Sequences) %>%
        dplyr::ungroup(.)

    abund_df <- primer_draws %>%
        select(Sequences, Primer_group) %>%
        mutate(Sequences = strsplit(Sequences, ",")) %>%
        unnest(Sequences)

    abund_matrix <- table(abund_df$Primer_group, abund_df$Sequences)

    out_matrix <- abund_matrix != 0

    if (nrow(out_matrix) > 1 && ncol(out_matrix) > 1) {
        out_matrix <- out_matrix[, colSums(out_matrix) > 0]
        out_matrix <- out_matrix[rowSums(out_matrix) > 0, ]
    } else if (nrow(out_matrix) > 1 && ncol(out_matrix) <= 1) {
        out_matrix <- out_matrix[rowSums(out_matrix) > 0, ]
    } else if (nrow(out_matrix) <= 1 && ncol(out_matrix) > 1) {
        out_matrix <- out_matrix[, colSums(out_matrix) > 0]
    }

    excluded_seqs <- filter(ag_data[[1]], !(Id %in% colnames(out_matrix)))

    Primer_candidates <- filter(primer_draws, primer_draws$Primer_group %in% rownames(out_matrix))

    message("Done!\n")

    new_prider(list(Description = description, Conversion_table = ag_data[[1]], Primer_candidates = primer_draws,
        Excluded_sequences = excluded_seqs, Primer_matrix = out_matrix))
}

#' @title print.prider
#' @rdname prider
#' @param x An object from prider function.
#' @param ... Other arguments.
#' @export
#' @importFrom dplyr count
print.prider <- function(x, ...) {
    descr <- x$Description
    input_seqs <- nrow(x$Conversion)
    excl_seqs <- nrow(x$Excluded_sequences)
    incl_seqs <- input_seqs - excl_seqs
    primer_candidates <- sum(x$Primer_candidates$Primer_group_size)
    groups <- nrow(dplyr::count(x$Primer_candidates, Primer_group))
    total <- paste(input_seqs, "input sequences;", incl_seqs, "included and", excl_seqs,
        "excluded.\n")
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
#' @param x An object from prider function.
#' @param ... Other arguments.
#' 
#' @export
#' @importFrom gplots heatmap.2
plot.prider <- function(x, ...) {
    matr <- x$Primer_matrix * 1
    if (ncol(matr) >= 2 && nrow(matr) >= 2) {
        gplots::heatmap.2(matr, scale = "none", trace = "none", col = c("white",
            "black"), xlab = "Sequence Id", ylab = "Primer cluster", key = FALSE)
    } else {
        stop("Primer_matrix too small to be plotted.")
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
#' @param prider_obj An object from prider function.
#' @param ix A number. The number of the primer cluster.
#' @return primer_obj
#' @examples
#'
#' test_fasta <- system.file('extdata', 'test.fasta', package = 'prider')
#'
#' primer_designs <- prider(test_fasta)
#'
#' primers(primer_designs)
#'
#' primers(primer_designs)[1]
#'
#' @export
primers <- function(prider_obj) UseMethod("primers")


#' @rdname primers
#' @export
primers.default <- function(prider_obj) warning(paste("Function 'primer' does not know how to handle object of class",
    class(prider_obj), "and can only be used on class 'prider'."))


#' @rdname primers
#' @export
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr ungroup
primers.prider <- function(prider_obj) {
    primer_obj <- prider_obj$Primer_candidates %>%
        dplyr::select(Primer_group, Seq_group_size, Primer_group_size, Primers) %>%
        dplyr::group_by(Primer_group) %>%
        dplyr::mutate(Primers = stringr::str_split(Primers, ",")) %>%
        dplyr::ungroup(.)
    new_primers(primer_obj)
}


#' @rdname primers
#' @param x An object from sequence function.
#' @param ... Other arguments.
#' @export
#' @importFrom dplyr select
print.primers <- function(x, ...) {
    cat("To access the primers within a group,\n")
    cat("please use subsetting, eg. primers(.)[42]\n")
    cat("\n")
    class(x) <- "data.frame"
    x %>%
        dplyr::select(Primer_group, Seq_group_size, Primer_group_size) %>%
        data.frame %>%
        print(.)
}


#' @rdname primers
#' @param primer_obj An object from sequence function.
#' @param ix A number. The number of the primer cluster.
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
#' @param prider_obj An object from prider function.
#' @param ix A number. The number of the primer cluster.
#' @return sequence_obj
#' @examples
#'
#' test_fasta <- system.file('extdata', 'test.fasta', package = 'prider')
#'
#' primer_designs <- prider(test_fasta)
#'
#' sequences(primer_designs)
#'
#' sequences(primer_designs)[1]
#'
#' @export
sequences <- function(prider_obj) UseMethod("sequences")


#' @rdname sequences
#' @export
sequences.default <- function(prider_obj) warning(paste("Function 'sequences' does not know how to handle object of class",
    class(x), "and can only be used on class 'prider'."))


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
    conversion_tbl <- prider_obj$Conversion_table
    sequence_obj <- prider_obj$Primer_candidates %>%
        dplyr::select(Primer_group, Seq_group_size, Primer_group_size, Sequences) %>%
        dplyr::group_by(Primer_group) %>%
        dplyr::mutate(Sequences = stringr::str_split(Sequences, ",")) %>%
        dplyr::ungroup(.) %>%
        tidyr::unnest(Sequences) %>%
        dplyr::left_join(conversion_tbl, by = c(Sequences = "Id")) %>%
        dplyr::select(-Sequences) %>%
        tidyr::nest(Sequences = c(Original_id, Sequence))
    new_sequences(sequence_obj)
}


#' @rdname sequences
#' @param x An object from sequence function.
#' @param ... Other arguments.
#' @export
#' @importFrom dplyr select
print.sequences <- function(x, ...) {
    cat("To access the sequences within a group,\n")
    cat("please use subsetting, eg. sequences(.)[42]\n")
    cat("\n")
    class(x) <- "data.frame"
    x %>%
        dplyr::select(Primer_group, Seq_group_size, Primer_group_size) %>%
        data.frame %>%
        print(.)
}


#' @rdname sequences
#' @param sequence_obj An object from sequence function.
#' @param ix A number. The number of the primer cluster.
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

