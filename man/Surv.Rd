% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{Surv}
\alias{Surv}
\title{Survival Function Wrapper}
\usage{
Surv(
  time,
  time2,
  event,
  type = c("right", "left", "interval", "counting", "interval2", "mstate"),
  origin = 0
)
}
\arguments{
\item{time}{
    for right censored data, this is the follow up time.  For interval
    data, the first argument is the starting time for the interval. 
  }

\item{time2}{
    ending time of the interval for interval censored  or counting
    process data only.  Intervals are assumed to be open on the left and
    closed on the right, \code{(start, end]}.  For counting process
    data, \code{event} indicates whether an event occurred at the end of
    the interval.
  }

\item{event}{
    The status indicator, normally 0=alive, 1=dead.  Other choices are
    \code{TRUE}/\code{FALSE} (\code{TRUE} = death) or 1/2 (2=death). For
    interval censored data, the status indicator is 0=right censored,
    1=event at \code{time}, 2=left censored, 3=interval censored.
    For multiple endpoint data the event variable will be a factor,
    whose first level is treated as censoring.
    Although unusual, the event indicator can be omitted, in which case
    all subjects are assumed to have an event.
  }

\item{type}{
    character string specifying the type of censoring. Possible values
    are \code{"right"}, \code{"left"}, \code{"counting"},
    \code{"interval"}, \code{"interval2"} or \code{"mstate"}.  
  }

\item{origin}{
    for counting process data, the hazard function origin.  This option
    was intended to be used in conjunction with a model containing
    time dependent
    strata in order to align the subjects properly when they cross over
    from one strata to another, but it has rarely proven useful.}
}
\description{
This alias provides direct access to the survival::Surv function.
}
