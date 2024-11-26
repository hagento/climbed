% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processArgs.R
\name{processArgs}
\alias{processArgs}
\title{Process and Standardize Input Arguments}
\usage{
processArgs(tLim, std, ssp)
}
\arguments{
\item{tLim}{A numeric value, numeric vector, or named list/vector with names
\code{"HDD"} and/or \code{"CDD"}.
 \itemize{
   \item If a single numeric value is provided, sequences for both \code{"HDD"}
   and \code{"CDD"} are generated by varying the value by -5 and +5.
   \item If a named list/vector is provided, each element is processed:
     \itemize{
       \item If an element is a single numeric value, it is converted into a
       sequence by varying it by -5 and +5.
       \item If an element is already a numeric sequence, it remains unchanged.
     }
 }}

\item{std}{A numeric value or a named numeric vector with names \code{"tLim"} and \code{"tAmb"}.
\itemize{
  \item If a single numeric value is provided, it is applied to both \code{"tLim"} and \code{"tAmb"}.
  \item If a named numeric vector is provided, it must have names \code{"tLim"} and \code{"tAmb"}.
}}

\item{ssp}{A character string or vector specifying the SSP (Shared Socioeconomic Pathways) scenarios.}
}
\value{
A list containing standardized arguments:
  \describe{
    \item{\code{tLim}}{A named list with \code{"HDD"} and \code{"CDD"} components, each being numeric sequences.}
    \item{\code{std}}{A named numeric vector with names \code{"tLim"} and \code{"tAmb"}.}
    \item{\code{ssp}}{A character vector of SSP scenarios.}
  }
}
\description{
The function processes and standardizes the input arguments \code{tLim}, \code{std},
and \code{ssp}, ensuring they are in the required format for further processing in
the \code{start} function or other applications.
}
\author{
Hagen Tockhorn
}
