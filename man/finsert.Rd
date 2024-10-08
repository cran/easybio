% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc.R
\name{finsert}
\alias{finsert}
\title{Insert Specific Values into a Character Vector at Defined Positions}
\usage{
finsert(
  x = expression(c(0, 1, 3) == "Neutrophil", c(2, 4, 8) == "Macrophage"),
  len = integer(),
  setname = TRUE,
  na = "Unknown"
)
}
\arguments{
\item{x}{An expression defining the value to insert and the positions at which
to insert them. The expression should be a list of logical comparisons,
where the left side is a numeric vector of positions and the right side is
the corresponding character value to insert.}

\item{len}{The desired length of the output character vector. If the specified
positions exceed this length, the vector will be padded with the \code{na} value.}

\item{setname}{A logical value indicating whether to set names for the elements
of the vector. If \code{TRUE}, names are set as character representations of the
positions from 0 to the length of the vector minus one.}

\item{na}{The default value to use for positions not specified in \code{x}. This
value is also used to pad the vector if its length exceeds the positions
specified in \code{x}.}
}
\value{
A named character vector with the specified values inserted at given
positions and padded with the \code{na} value if necessary.
}
\description{
This function constructs a character vector of a specified length, inserting
given values at positions determined by numeric indices. It is designed for
single cell annotation tasks, where specific annotations need to be placed
at certain positions in a vector.
}
\examples{
# Example usage:
# Insert "Neutrophil" at positions 0, 1, 3 and "Macrophage" at positions 2, 4, 8
# in a vector of length 10, with "Unknown" as the default value.
library(easybio)
annotated_vector <- finsert(
  x = expression(
    c(0, 1, 3) == "Neutrophil",
    c(2, 4, 8) == "Macrophage"
  ),
  len = 10,
  na = "Unknown"
)
print(annotated_vector)
}
