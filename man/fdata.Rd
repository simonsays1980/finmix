\name{fdata}
\docType{class}

\alias{class:fdata}
\alias{fdata}

%accessor
\alias{getY, fdata-method}
\alias{setY<-, fdata-method}
\alias{getColY, fdata-method}
\alias{getRowY, fdata-method}
\alias{getN, fdata-method}
\alias{setN<-, fdata-method}
\alias{getR, fdata-method}
\alias{setR<-, fdata-method}
\alias{getS, fdata-method}
\alias{setS<-, fdata-method}
\alias{getColS, fdata-method}
\alias{getRowS, fdata-method}
\alias{getBycolumn, fdata-method}
\alias{setBycolumn<-, fdata-method}
\alias{getName, fdata-method}
\alias{setName<-, fdata-method}
\alias{getType, fdata-method}
\alias{setType<-, fdata-method}
\alias{getSim, fdata-method}
\alias{setSim<-, fdata-method}
\alias{getExp, fdata-method}
\alias{setExp<-, fdata-method}
\alias{getColExp, fdata-method}
\alias{getRowExp, fdata-method}
\alias{getY, fdata-method}
\alias{setY<-, fdata-method}
\alias{getT, fdata-method}
\alias{setT<-, fdata-method}
\alias{getColT, fdata-method}
\alias{getRowT, fdata-method}

% constructor
\alias{fdata}

% checking
\alias{hasY, fdata-method}
\alias{hasS, fdata-method}
\alias{hasExp, fdata-method}
\alias{hasT, fdata-method}

% show
\alias{show, fdata-method}

% plot
\alias{plot, fdata-method}

\title{Finmix Data}
\description{
    The \code{fdata} class stores the data for finite mixture distributions.
}

\details{
    The \code{fdata} class is constructed by calling its constructor 
    \code{fdata()}. All arguments in the constructor are optional. 
}

\section{Constructor}{
    \describe{\code{fdata(y = matrix(), N = 1, r = 1, S = matrix(),
            bycolumn = TRUE, name = character(), type = "discrete", 
            sim = FALSE, exp = matrix(), T = matrix())}:

        Constructs an \code{fdata} object from the input arguments. All
        arguments are optional. If provided, arguments \code{y}, \code{S}, 
        \code{exp} and \code{T} must be of type \code{matrix}. 

        To construct an empty \code{fdata} object the constructor can be
        without any argument provided: \code{fdata()}. 
    }
}

\section{Accessors}{
In the following code snippets, \code{x} is an \code{fdata} object and the symbol
\code{@} represents a slot of this \code{S4} object.
    \describe{
        \item{}{\code{getY(x)}, \code{setY(x)<-}:
            Get and set the \code{matrix} containing the observations 
            \code{@y} in \code{x}. In case the setter is called the number
            of observations \code{@N} and the variable dimension \code{r}
            are set automatically.  
        }
        \item{}{\code{getColY(x)}, \code{getRowY(x)}:
            Get the \code{matrix} of observations in \code{@y} of \code{x} either
            sorted by row or by column.            
        }
        \item{}{\code{getS(x)}, \code{setS(x)<-}:
            Get and set the indicator \code{matrix} in \code{@S} of \code{x}. 
            Indicators are stored as \code{integer} and must be of dimension 
            \code{@N x 1} in case of \code{@bycolumn = TRUE} or of dimension
            \code{1 x @N} in case of \code{@bycolumn = FALSE}. Certain safe
            guards check for consistency with the remaining slots of \code{x}.
        }
        \item{}{\code{getColS(x)}, \code{getRowS(x)}:
            Get the \code{matrix} of indicators in \code{@S} of \code{x} either 
            sorted by row or by column.
        }
        \item{}{\code{getN(x)}, \code{setN(x)<-}:
            Get and set the number of observations in \code{@N} of the
            \code{matrix} of observations in \code{@y}. The number of 
            observations is stored as an \code{integer} and certain safe guards
            check for consistency with the remaining slots of \code{x}.    
        }
        \item{}{\code{getR(x)}, \code{setR(x)<-}:
            Get and set the dimension of variables in \code{@r} of 
            the \code{matrix} of observations in \code{@y} of \code{x}.
            The dimension of variables is stored as an \code{integer} and certain
            safe guards check for consistency with the remaining slots of 
            \code{x}. For univariate distributions \code{@r} must be \code{1}
            and for multivariate distributions \code{@r} must be \code{>1}.         
        }
        \item{}{\code{getBycolumn(x)}, \code{setBycolumn(x)<-}:
            Get and set the ordering \code{@bycolumn} of observations in the 
            \code{matrix} of \code{@y} in \code{x}. The ordering is stored
            as a \code{logical} and must be either \code{TRUE} or \code{FALSE}.
            If \code{@bycolumn} is set to \code{TRUE} observations in the 
            \code{matrix} of \code{@y} are ordered in columns, otherwise 
            observations are ordered in rows. Certain safe guards check for
            consistency with the remaining slots.
        }
        \item{}{\code{getType(x)}, \code{setType(x)<-}:
            Get and set the data type \code{@type} of the observations in the
            \code{matrix} of \code{@y}. The data type is stored as \code{character}
            and it must be either \code{discrete} or \code{continuous}. The
            data type determines the behavior of certain function calls.          
        }
        \item{}{\code{getSim(x)}, \code{setSim(x)<-}:
            Get and set the simulation indicator in \code{@sim}. This slot 
            indicates if the observations in the \code{matrix} of \code{@y}
            were simulated. The simulation indicator is stored as \code{logical}
            and must be either \code{TRUE} or \code{FALSE}.        
        }
        \item{}{\code{getName(x)}, \code{setName(x)<-}:
            Get and set the name of the data in \code{@name} of \code{x}. The
            name is optional and set as \code{character()} if no name is
            provided by the user.       
        }
        \item{}{\code{getExp(x)}, \code{setExp(x)<-}:
            Get and set the \code{matrix} of data exposures \code{@exp} in \code{x}. 
            This is optional and becomes only relevant in case a \code{poisson} 
            mixture is fitted to the data. Certain safe guards check for
            consistency with the remaining slots of \code{x}. In case \code{@y},
            \code{@S} and \code{@T} are empty, the number of observations \code{@N}
            and the dimension of variables \code{@r} are set automatically.            
        }
        \item{}{\code{getColExp(x)}, \code{getRowExp(x)}:
            Get the \code{matrix} of exposures in \code{@exp} of \code{x} either
            sorted by row or by column.
        }
        \item{}{\code{getT(x)}, \code{setT(x)<-}:
            Get and set the repetitions \code{matrix} in \code{@T} of \code{x}. 
            Repetitions are optional and become only relevant in case a 
            \code{binomial} mixture is fitted to the data. Repetitions are stored as
            \code{integer} in \code{x}. Certain safe guards check
            for consistency with the remaining slots of \code{x}. In case \code{@y},
            \code{@S} and \code{@exp} are emoty, the number of observations \code{@N}
            and the dimension of variables \code{@r} are set automatically.
        }
        \item{}{\code{getColT(x)}, \code{getRowT(x)}:
            Get the \code{matrix} of repetitions in \code{@T} of \code{x} either
            sorted by row or by column.
        }
    }
}

\section{Checking}{
    In the following code snippets, \code{x} is an \code{fdata} object and the symbol
    \code{@} represents a slot of this \code{S4} object.
    \describe{
        \item{}{\code{hasY(x, verbose = FALSE)}:
            Checks wether \code{@y} of \code{x} is empty. Returns \code{TRUE}, if
            \code{@y} contains \code{matrix} with not all entries \code{NA}, 
            otherwise it returns \code{FALSE}. If \code{verbose} is set to 
            \code{TRUE} an error is thrown in case of \code{@y} being empty. 
        }
        \item{}{\code{hasS(x, verbose = FALSE)}:
            Checks wether \code{@S} of \code{x} is empty. Returns \code{TRUE}, if
            \code{@S} contains \code{matrix} with not all entries \code{NA}, 
            otherwise it returns \code{FALSE}. If \code{verbose} is set to 
            \code{TRUE} an error is thrown in case of \code{@S} being empty. 
        }
        \item{}{\code{hasExp(x, verbose = FALSE)}:
            Checks wether \code{@exp} of \code{x} is empty. Returns \code{TRUE}, if
            \code{@exp} contains \code{matrix} with not all entries \code{NA}, 
            otherwise it returns \code{FALSE}. If \code{verbose} is set to 
            \code{TRUE} an error is thrown in case of \code{@exp} being empty. 
        }
        \item{}{\code{hasT(x, verbose = FALSE)}:
            Checks wether \code{@T} of \code{x} is empty. Returns \code{TRUE}, if
            \code{@T} contains \code{matrix} with not all entries \code{NA}, 
            otherwise it returns \code{FALSE}. If \code{verbose} is set to 
            \code{TRUE} an error is thrown in case of \code{@T} being empty. 
        }
    }
}

\section{Plotting}{
    \describe{\code{plot(x, dev = TRUE)}:

        Plots the observations in \code{@y} of \code{x}. In case the data 
        type \code{@type} of \code{x} is \code{discrete} a \code{\link{barplot}} 
        is used for visualizing. In case the data type is \code{continuous} a
        histogram is plotted by using \code{\link{hist}}. If the dimension of 
        variables \code{@r} is equal to \code{2} a bivariate Kernel density is
        plotted using \code{\link{bkde2D}} from package \code{\link{KernSmooth}}
        for the Kernel density estimates, \code{\link{persp}} for a persepctive
        plot of the denisty and \code{\link{contour}} for a corresponding 
        contour plot. In case the dimension of variables \code{@r} is \code{>2},
        histograms are plotted for all variables and \code{\link{pairs}} is used
        for scatterplots of all pairs of data.

        If argument \code{dev = FALSE} no graphical device is opened and the 
        user is able to store all plots to a file using \code{\link{pdf}},
        \code{\link{png}}, etc.   
    }
}
\author{ Lars Simon Zehnder }
\examples{
    fdata.obj <- fdata(y = rpois(100, 312), sim = TRUE)
    fdata.obj
}
\keyword{classes}
\keyword{methods}
