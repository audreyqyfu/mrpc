\name{Cut_Modules}
\alias{Cut_Modules}

\title{
Cut the Modules
}
\description{
Similar to \link{cut2} function with some modification.
}
\usage{
Cut_Modules(x, cuts, m, g, levels.mean = FALSE,digits,minmax=TRUE,
oneval=TRUE, onlycuts=FALSE)

}

\arguments{
  \item{x}{
numeric vector to classify into intervals
}
  \item{cuts}{
cut points
}
  \item{m}{
desired minimum number of observations in a group
}
\item{g}{
number of quantile groups
}
\item{levels.mean}{
set to TRUE to make the new categorical vector have levels attribute that is the group means of x instead of interval endpoint labels
}
\item{digits}{
number of significant digits to use in constructing levels. Default is 3 (5 if levels.mean=TRUE)
}

\item{minmax}{
if cuts is specified but min(x)<min(cuts) or max(x)>max(cuts), augments cuts to include min and max x
}
\item{oneval}{
if an interval contains only one unique value, the interval will be labeled with the formatted version of that value instead of the interval endpoints, unless oneval=FALSE
}
\item{onlycuts}{
set to TRUE to only return the vector of computed cuts. This consists of the interior values plus outer ranges.
}

}

\value{
Vector
}
\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

