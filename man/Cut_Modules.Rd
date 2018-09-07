\name{Cut_Modules}
\alias{Cut_Modules}

\title{
Cut the Modules
}
\description{
Similar to \link{cut2} function with some modification.
}
\usage{
Cut_Modules(x, cuts, m, g, levels.mean = FALSE, digits, minmax = TRUE,
    oneval = TRUE, onlycuts = FALSE)

}

\arguments{
  \item{x}{
Numeric vector to classify into intervals.
}
  \item{cuts}{
Cut points
}
  \item{m}{
Desired minimum number of observations in a group.
}
\item{g}{
Number of quantile groups.
}
\item{levels.mean}{
Set to TRUE to make the new categorical vector have levels attribute that is the group means of x instead of interval endpoint labels.
}
\item{digits}{
Number of significant digits to use in constructing levels. The default is 3, and 5 if levels.mean = TRUE.
}

\item{minmax}{
If cuts is specified but min(x) < min(cuts) or max(x) > max(cuts) augments cuts to include min and max x.
}
\item{oneval}{
If an interval contains only one unique value, the interval will be labeled with the formatted version of that value instead of the interval endpoints unless oneval = FALSE.
}
\item{onlycuts}{
Set to TRUE to only return the vector of computed cuts. This consists of the interior values plus outer ranges.
}

}

\value{
Vector
}
\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>, Evan A Martin <mart9986@vandals.uidaho.edu>, and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}

