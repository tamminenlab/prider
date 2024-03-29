% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{chunker}
\alias{chunker}
\title{chunker}
\usage{
chunker(seq_table, window_size = 20L)
}
\arguments{
\item{seq_table}{A DataFrame containing a column for sequence ids (Id) and
sequences (Seq).}

\item{window_size}{An integer. Set the sliding window width.}
}
\value{
A DataFrame containing columns for the sequence ids (Id),
        indexes (Ix), joined ids and indexes (Id_Ix), and the primer
        sequences (Seq).
}
\description{
Creates all primer candidates for a group of sequences using a
             sliding window.
}
\details{
Sliding window to create chunks of DNA sequences
}
\examples{

test_csv <- system.file("extdata", "test.csv", package = "prider")

test_csv <- read.csv(test_csv)

chunks <- chunker(test_csv)

}
