\name{model}
\docType{class}

\alias{class:model}
\alias{model}

%accessor
\alias{getDist, model-method}
\alias{setDist<-, model-method}
\alias{getR, model-method}
\alias{setR<-, model-method}
\alias{getWeight, model-method}
\alias{setWeight<-, model-method}
\alias{getPar, model-method}
\alias{setPar<-, model-method}
\alias{getIndicmod, model-method}
\alias{setIndicmod<-, model-method}
\alias{getIndicfix, model-method}
\alias{setIndicfix<-, model-method}
\alias{getT, model-method}
\alias{setT<-, model-method}

% constructor
\alias{model}

% checking
\alias{hasWeight, model-method}
\alias{hasPar, model-method}
\alias{hasT, model-method}

% show
\alias{show, model-method}

% plot
\alias{plot, model-method}
\alias{plotPointProc, model-method}

% tools
\alias{mixturemar, model-method}

\title{Finmix Model}
\description{
    The \code{model} class the model for a finite mixture distribution.
}

\details{
    The \code{model} class is constructed by calling its constructor 
    \code{model()}. All arguments in the constructor are optional. 
}

\section{Constructor}{
    \describe{\code{model(dist = "poisson", r, K, weight = matrix(),
            par = list(), indicmod = "multinomial", indicfix = TRUE,
            T = matrix())}:

        Constructs a \code{model} object from the input arguments. All
        arguments are optional. If provided, argument \code{T} must be 
        of type \code{matrix}.  

        To construct an empty \code{model} object the constructor can be
        without any argument provided: \code{model()}. 
    }
}

\section{Accessors}{
In the following code snippets, \code{x} is an \code{model} object and the symbol
\code{@} represents a slot of this \code{S4} object.
    \describe{
        \item{}{\code{getDist(x)}, \code{setDist(x)<-}:
            Get and set the distribution of the mixture model in \code{@dist} 
            of \code{x}. The following
            distributions are implemented: \code{"poisson"}, \code{"binomial"},
            \code{"exponential"}, \code{"normal"}, \code{"student"}, 
            \code{"normult"} (multivariate Normal) and \code{"studmult"} 
            (multivariate Student-t). Only models with the same distributions
            are implemented. 
        }
        \item{}{\code{getR(x)}, \code{setR(x)<-}:
            Get and set the dimension of variables in \code{@r} of \code{x}. 
            The dimension of variables is stored as an \code{integer} and certain
            safe guards check for consistency with the remaining slots of 
            \code{x}. For univariate distributions (\code{"poisson"}, 
            \code{"binomial"}, \code{"exponential"}, \code{"normal"}, 
            \code{"student"}) \code{@r} must be \code{1} and for multivariate 
            distributions (\code{"normult"}, \code{"studmult"}) \code{@r} 
            must be \code{>1}.         
        }
        \item{}{\code{getK(x)}, \code{setK(x)<-}:
            Get and set the number of components in \code{@K} of \code{x}.
            The number of components must be a positive integer.
            It is stored as an \code{integer} and certain safe guards check 
            for validity. 
        }
        \item{}{\code{getWeight(x)}, \code{setWeight(x)<-}:
            Get and set the weights of the finite mixture model in \code{@weight}
            of \code{x}. The weights must be a \code{1 x @K} \code{matrix} of type 
            \code{numeric} all \code{<1} and \code{>0} and must sum to \code{1}. 
            Certain safe guards check for validity and consistency with the remaining 
            slots of \code{x}.
        }
        \item{}{\code{getPar(x)}, \code{setPar(x)<-}:
            Get and set the component parameters of the finite mixture model
            in \code{@par} of \code{x}. If the setter is called parameters 
            must be provided in a \code{list} with appropriate naming regarding
            the distribution of the model in slot \code{@dist}. The following 
            naming rules apply:
            \itemize{
                \item \code{"poisson"}: A \code{vector} of positive Poisson 
                    parameters with name \code{$lambda} in the \code{list}.           
                \item \code{"binomial"}: A \code{vector} of positive Binomial
                    parameters with name \code{$p} in the \code{list}. All 
                    parameters must be \code{>0} and \code{<1}.
                \item \code{"exponential"}: A \code{vector} of positive 
                    Exponential parameters with name \code{$lambda} in the 
                    \code{list}.
                \item \code{"normal"}: A \code{vector} of means with name
                    \code{$mu} in the \code{list} and a \code{vector} of 
                    standard deviations with name \code{$sigma} in the \code{list}. 
                    All standard deviations must be positive.
                \item \code{"student"}: A \code{vector} of location parameters
                    with name \code{$mu} in the \code{list} and a \code{vector} 
                    of scale parameters with name \code{$sigma} in the \code{list}. 
                    All scale parameters must be positive. In addition the 
                    degrees of freedom must be provided as a \code{vector}
                    with name \code{$df} in the \code{list}.
                \item \code{"normult"}: An \code{array} or \code{matrix} of 
                    dimension \code{@r x @K} containing the means for each dimension
                    and component named \code{$mu} in the \code{list}.
                    Further, an \code{array} of dimension \code{@r x @r x @K}
                    containing the variance-covariance matrices named \code{$sigma}
                    in the \code{list}. All matrices must be stored as a \code{matrix}
                    and must be positive-definite. 
                \item \code{"studmult"}: An \code{array} or \code{matrix} of 
                    dimension \code{@r x @K} containing the location parameters
                    for each dimension and component named \code{$mu} in the \code{list}.
                    Further, an \code{array} of dimension \code{@r x @r x @K}
                    containing the scale matrices named \code{$sigma}
                    in the \code{list}. All matrices must be stored as a \code{matrix}
                    and must be positive-definite. In addition, degrees of freedom
                    must be provided as a \code{vector} with name \code{$df} in
                    the \code{list}.                   
            }
        }
        \item{}{\code{getIndicmod(x)}, \code{setIndicmod(x)<-}:
            Get and set the indicator model in \code{@indicmod} of \code{x}. 
            Each finite mixture model has an underlying model for its indicators.
            Right now only the model \code{"multinomial"} is implemented.       
        }
        \item{}{\code{getIndicfix(x)}, \code{setIndicfix(x)<-}:
            Get and set the indicator for a model with fixed indicators in 
            \code{@indicfix} of \code{x}. A finite mixture model can have 
            predefined indicators, either because they are observed or 
            estimated by pre-sample classification. This indicator slot
            is stored as \code{logical} and must be either \code{TRUE} in
            case fixed indicators are provided in the \code{\link{fdata}} 
            or \code{FALSE} if otherwise.
        }
        \item{}{\code{getT(x)}, \code{setT(x)<-}:
            Get and set the repetitions \code{matrix} in \code{@T} of \code{x}. 
            Repetitions are optional and become only relevant in case the 
            distribution in \code{@dist} is set to \code{"binomial"}. 
            Repetitions must be stored in a \code{matrix} in case the 
            setter is called.            
        }
    }
}

\section{Checking}{
    In the following code snippets, \code{x} is an \code{model} object and the symbol
    \code{@} represents a slot of this \code{S4} object.
    \describe{
        \item{}{\code{hasWeight(x, verbose = FALSE)}:
            Checks wether \code{@weight} of \code{x} is empty. Returns \code{TRUE}, 
            if \code{@weight} contains a \code{matrix} with not all entries 
            \code{NA} and dimension \code{1 x @K}, otherwise it returns \code{FALSE}.
            If \code{verbose} is set to \code{TRUE} an error is thrown in case of 
            \code{@weight} being empty or having a wrong dimension.
        }
        \item{}{\code{hasPar(x, verbose = FALSE)}:
            Checks wether \code{@par} of \code{x} is empty. Returns \code{TRUE}, if
            \code{@par} contains a \code{list} with length appropriately named 
            entries of correct dimension. See \code{\link{setPar}} for defining
            correct parameters for a finite mixture model. In case \code{@par}
            is either empty or incorrectly specified the function returns \code{FALSE}.  
            If \code{verbose} is set to \code{TRUE} an error is thrown in case of 
            \code{@par} being empty or incorrectly specified. 
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
        
        Plots a model as specified by \code{x}. 
        The following types of plots are returned in regard to \code{@dist}:
        \itemize{
            \item \code{"poisson"}: A barplot of the probabilities 
                over the range between minimum and maximum of the Poisson mixture 
                distribution.
            \item \code{"binomial"}: A line diagram of the probabilities over 
                the range between and minimum and maximum of the Binomial mixture 
                distribution.
            \item \code{"exponential"}: An density plot over the range 
                between minimum and maximum of the Exponential mixture 
                distribution.
            \item \code{"normal"}: A density plot over the range between minimum
                and maximum of the Normal mixture distribution. 
            \item \code{"student"}: A density plot over the range between minimum
                and maximum of the Student-t distribution.
            \item \code{"normult"}: In case the dimension of variables in 
                \code{@r} is equal to \code{2} a perspective plot and a contour 
                plot of the bivariate density of the bivariate Normal mixture
                is returned. In case \code{@r} is \code{>2} contour plots for
                all possible marginal Normal mixture models of dimension \code{2} 
                are returned.
            \item \code{"studmult"}: In case the dimension of variables in 
                \code{@r} is equal to \code{2} a perspective plot and a contour 
                plot of the bivariate density of the bivariate Student-t mixture
                is returned. In case \code{@r} is \code{>2} contour plot for
                all possible marginal Student-t mixture models of dimension \code{2} 
                are returned.
        }
        If argument \code{dev = FALSE} no graphical device is opened and the 
        user is able to store all plots to a file using \code{\link{pdf}},
        \code{\link{png}}, etc.   
    }
    \describe{\code{plotPointProc(x, dev = TRUE, ...)}: 
        
        Plots the point process representation for the mixture model specified by 
        \code{x}. The following type of plots are returned in regard to \code{@dist}:
        \itemize{
            \item \code{"poisson"}: A point plot, indicating the position of the
                components and their corresponding weights by points of certain
                sizes.
        }
    }
}
\author{ Lars Simon Zehnder }
\examples{
    model.obj <- model(dist = "binomial", K = 2, indicfix = TRUE)
    model.obj
    setT(model.obj)         <- as.matrix(100)
    setPar(model.obj)       <- list(p = c(.3, .7))
    setWeight(model.obj)    <- matrix(c(.1, .9), nrow = 1, ncol = 2)
    plot(model.obj) 
}
\keyword{classes}
\keyword{methods}
