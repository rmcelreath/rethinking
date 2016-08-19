\name{Trolley}
\alias{Trolley}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Experimental data on ethical dilemmas}
\description{
  These data comprise outcomes of experimental ethical dilemmas that are often called 'trolley' problems. Data kindly provided by Fiery Cushman.
}
\usage{
data(Trolley)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item case: a code that combines treatment and story labels
        \item response: participant's rating of appropriateness of action in story
        \item order: integer order story appeared, within participant
        \item id: participant id (factor)
        \item age: participant's age in years
        \item male: participant's gender; 1 for male, 0 for female
        \item edu: participant's highest educational level
        \item action: treatment code for story with action (1) or inaction (0)
        \item intention: treatment code for intent (1) or lack of intent (0)
        \item contact: treatmetn code for contact action (1) or lack of contact (0)
        \item story: factor label for basic scenario modified by treatments
        \item action2: alternative coding of action that is union of action and contact variables
    }
}
\value{
}
\references{Cushman et al. 2006. Psychological Science 17:1082--1089.}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

