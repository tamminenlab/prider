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


#' Checks if set1 is a subset of set2. The sets are comma-delimited strings.
#'
#' @param set1 A string
#' @param set2 A string
#' @return Boolean
is_subset_of <- function(set1, set2)
{
    set1 <- str_split(set1, ",") %>% unlist
    set2 <- str_split(set2, ",") %>% unlist
    all(set1 %in% set2)
}


#' Runs the Prider workflow
#'
#' @param input_fasta A DataFrame or a string
#' @param primer_length A number
#' @return A list containing sequence names, sequence network, primer clusters,
#'         sequence clusters, cluster overlap and primer list
#' @examples
#'
#' primer_designs <- prider("test.fasta", primer_length = 20)
#' 
#' @export
#' @importFrom stringr str_length
#' @importFrom stringr str_split
#' @importFrom purrr pluck
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom purrr map2
#' @importFrom purrr map2_lgl
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components
#' @importFrom igraph groups
#' @importFrom tidyr separate
#' @importFrom tidyr unnest
#' @importFrom tidyr unite
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr as_tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom dplyr ungroup
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom dplyr row_number
#' @importFrom dplyr n
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom blaster blast
#' @importFrom blaster read_fasta
prider <- function(input_fasta,
            primer_length,
            min_identity = 0.95)
{
    if (min_identity > 1 || min_identity < 0.6)
        stop("min_identity should be between 0.6 and 1")

    if (is.character(input_fasta))
        input_fasta <- read_fasta(input_fasta)

    message("Preparing BLAST data.")
    fasta_table <-
        input_fasta %>% 
        dplyr::mutate(Seq_no = paste0("P", row_number()))

    conversion_table <- 
        fasta_table %>%
        select(Seq_no, Id) %>% 
        dplyr::rename(Original_id = Id, Id = Seq_no) %>% 
        dplyr::as_tibble(.)

    database <- 
        fasta_table %>%
        dplyr::select(Seq_no, Seq) %>%
        dplyr::rename(Id = Seq_no) %>% 
        dplyr::as_tibble(.)

    primers <- 
        database %>%
        chunker(window_size = primer_length) %>% 
        dplyr::as_tibble(.)

    clust_primers <- 
        primers %>% 
        group_by(Seq) %>%
        summarise(Id = paste0(Id, collapse = ";"))

    accepts <-
        database %>% 
        dim %>% 
        purrr::pluck(1)
    
    message("Running BLAST.")
    tmp_output <- 
        blast(query = clust_primers,
              db = database,
              minIdentity = min_identity,
              maxAccepts = accepts, 
              alphabet = "nt", 
              strand = "plus",
              output_to_file = TRUE)
    on.exit(if (exists(tmp_output)) file.remove(tmp_output), add = TRUE)

    seq_len <-
        primers$Seq %>%
        purrr::pluck(1) %>%
        stringr::str_length(.) - 1

    message("")
    message("Processing BLAST output.")
    data_list <- 
        process_blast_table(tmp_output,
                            hit_len = seq_len) %>% 
        purrr::map(as_tibble) %>% 
        setNames(c("primer_network", "seq_network"))

    message("Clustering primers.")
    primer_clusters <- 
        data_list$primer_network %>% 
        igraph::graph_from_data_frame(directed=FALSE) %>% 
        igraph::components(.) %>%
        igraph::groups(.) %>%
        purrr::map_df(~dplyr::tibble(Id=.), .id='Primer_cluster') %>% 
        dplyr::left_join(primers, by = "Id") %>% 
        dplyr::group_by(Primer_cluster) %>%
        dplyr::mutate(Seq = make_degenerate_sequence(Seq, sequence_length = seq_len + 1)) %>% 
        dplyr::ungroup(.)
    
    message("Clustering sequence targets.")
    seq_clusters <- 
        primer_clusters %>% 
        tidyr::separate(Id, c("Seq_id", "Locus"), sep="_") %>% 
        dplyr::group_by(Seq) %>%
        dplyr::summarise(
                   Seq_n = n(),
                   Seq_series = paste(sort(Seq_id), collapse=","),
                   .groups = "drop") %>%
        dplyr::group_by(Seq_series, Seq_n) %>%
        dplyr::summarise(
                   Primer_n = n(),
                   Primer_series = paste(sort(Seq), collapse=","),
                   .groups = "drop") %>%
        dplyr::arrange(desc(Seq_n))

    message("Removing redundant sub-clusters.")
    subsets <- 
        seq_clusters %>%
        dplyr::pull(Seq_series) %>% 
        combn(2) %>%
        t %>% 
        as.data.frame %>% 
        dplyr::as_tibble(colnames = c("V1", "V2")) %>% 
        dplyr::mutate(V1 = as.character(V1),
                      V2 = as.character(V2)) %>% 
        dplyr::bind_rows(., rename(., V1 = V2, V2 = V1)) %>% 
        dplyr::mutate(
                   Is_subset = map2_lgl(V1, V2, ~is_subset_of(.x, .y))) %>% 
        dplyr::filter(Is_subset == TRUE) %>% 
        dplyr::mutate(
                   V1_n = map_dbl(V1, ~length(unlist(str_split(., ",")))),
                   V2_n = map_dbl(V2, ~length(unlist(str_split(., ",")))))

    redundant_groups <- 
        subsets %>%
        dplyr::group_by(V1) %>%
        dplyr::mutate(Max = max(V2_n)) %>%
        dplyr::filter(V2_n == Max,
                      V1_n != V2_n) %>%
        dplyr::select(V1)

    seq_clusters<- 
        seq_clusters %>%
        dplyr::ungroup(.) %>% 
        dplyr::filter(!(Seq_series %in% redundant_groups$V1)) %>% 
        dplyr::mutate(Cluster = row_number()) %>% 
        dplyr::select(Cluster, Seq_n, Primer_n, Seq_series, Primer_series)

    cluster_ids <-
        seq_clusters %>% 
        dplyr::select(Cluster, Seq_series)
    
    message("Detecting cluster overlaps.")
    cluster_overlap <- 
        seq_clusters %>%
        dplyr::pull(Seq_series) %>% 
        combn(2) %>%
        t %>% 
        as.data.frame %>% 
        dplyr::as_tibble(.) %>% 
        dplyr::mutate(
                   V1 = as.character(V1),
                   V2 = as.character(V2),
                   X1 = purrr::map(V1, ~unlist(str_split(., ","))),
                   X2 = purrr::map(V2, ~unlist(str_split(., ","))),
                   Int = purrr::map2(X1, X2, ~intersect(.x, .y)),
                   Empty = purrr::map_lgl(Int, ~identical(., character(0))),
                   Overlap = purrr::map_chr(Int, ~paste0(., collapse = " ")),
                   Len_C1 = purrr::map_dbl(X1, length),
                   Len_C2 = purrr::map_dbl(X2, length),
                   Len_Int = purrr::map_dbl(Int, length), 
                   Perc_of_cluster1 = round(Len_Int / Len_C1 * 100, 2), 
                   Perc_of_cluster2 = round(Len_Int / Len_C2 * 100, 2)) %>%
        dplyr::filter(!Empty) %>% 
        dplyr::left_join(cluster_ids, by = c("V1" = "Seq_series")) %>% 
        dplyr::rename(Cluster1 = Cluster) %>% 
        dplyr::left_join(cluster_ids, by = c("V2" = "Seq_series")) %>% 
        dplyr::rename(Cluster2 = Cluster) %>% 
        dplyr::select(Cluster1, Cluster2, Overlap, Perc_of_cluster1, Perc_of_cluster2)

    message("Creating primer list.")
    primer_list <- 
        seq_clusters %>%
        dplyr::select(Cluster, Primer_series, Seq_series) %>% 
        dplyr::mutate(
                   Seq_series = map(Seq_series, ~unlist(str_split(., ",")))) %>% 
        tidyr::unnest(cols = c(Seq_series)) %>%
        dplyr::mutate(
                   Primer_series = map(Primer_series, ~unlist(str_split(., ",")))) %>%
        tidyr::unnest(cols = c(Primer_series)) %>% 
        dplyr::left_join(primer_clusters,
                         by = c("Primer_series" = "Seq")) %>% 
        dplyr::select(-Seq_series) %>%
        dplyr::select(Cluster, Primer_cluster, Id, Primer_series) %>% 
        dplyr::rename(Primer_sequence = Primer_series) %>% 
        group_by(Cluster, Primer_cluster, Primer_sequence) %>%
        summarise(Ids = paste0(Id, collapse = ","))
    
    list(conversion_table = conversion_table, 
         seq_network = data_list$seq_network, 
         primer_clusters = primer_clusters,
         seq_clusters = seq_clusters,
         cluster_overlap = cluster_overlap,
         primer_list = primer_list)
}


#' Prepares a Gephi network file from data_list
#'
#' @param data_list A list
#' @param output_gephi_file A string
#' @export
make_gephi <- function(data_list, 
                output_gephi_file)
{
    annotation <- 
        data_list$seq_clusters %>%
        dplyr::select(Seq_series) %>%
        unique %>% 
        dplyr::mutate(Seq_series = map(Seq_series, ~unlist(str_split(., ","))),
                      Clust_no = row_number()) %>%
        tidyr::unnest(cols = c(Seq_series)) %>%
        tidyr::unite(Annot, Seq_series, Clust_no, sep=",") %>%
        dplyr::pull(Annot)

    seq_netw <-
        data_list$seq_network %>% 
        tidyr::unite(Netw, Left, Right, sep=",") %>%
        dplyr::pull(Netw)

    gdf <- c("nodedef>name VARCHAR,type VARCHAR",
             annotation,
             "edgedef>node1 VARCHAR,node2 VARCHAR",
             seq_netw)

    write(gdf, output_gephi_file)
}
