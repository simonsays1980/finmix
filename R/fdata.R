## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with finmix. If not, see <http://www.gnu.org/licenses/>.

#' Finmix `fdata` class
#' 
#' @description
#' The `fdata` class holds the data for finite mixture distributions.
#' 
#' @details 
#' The `fdata` class defines an essential part of the `finmix` package and 
#' MCMC sampling for finite mixture distributions. It stores the data for 
#' finite mixture distributions which includes always the observations stored 
#' in slot `y` and occasionally also known indicators in slot `S`. The latter 
#' ones define either a so-called finite mixture model with *fixed* indicators 
#' or are used as starting indicators in MCMC sampling for a model with unknown 
#' indicators. 
#' 
#' Observations can be stored in either in row or column format (default). In 
#' the former case the slot `bycolumn` has to be set to `FALSE` to indicate the 
#' safeguard functions in methods that the observations are stored in row 
#' format. If indicators are stored in the `fdata` object they must be stored 
#' in the same format as the observations. When using the setter `setS()<-` 
#' converting the repetitions to the right format is handled for the user. 
#' 
#' For discrete mixture models with Poisson or Exponential distributions 
#' exposures can be added to the data (and model). Exposures scale the rate 
#' parameters individually for each observation. Exposures get stored in the 
#' slot `exp` and have to be either of dimension `Nx1` or of dimension `1x1`. 
#' Like observations and indicators, exposures also have to be provided in the 
#' same data format, i.e. either row or column depending on the slot `bycolumn` 
#' set to `FALSE` or `TRUE`. When using the setter `setExp()<-` converting the 
#' repetitions to the right format is handled for the user. 
#' 
#' For mixtures of binomial distributions it is possible to include repetitions 
#' in the slot `T` of the `fdata` object. Repetitions can be constant or 
#' varying. In the former case the dimension of slot `T` is `1x1` and in the 
#' latter one it is `Nx1`. Depending on the slot `bycolumn` the repetitions 
#' have to be provided in row or column format. When using the setter 
#' `setT()<-` converting the repetitions to the right format is handled for the 
#' user. 
#' 
#' For mixtures of multivariate data the slot `r` is larger than one. For all 
#' other mixtures it is equal to one. Note that in case of multivariate mixture 
#' models the data in slot `y` has to be of dimension `Nxr` or `rxN` depending 
#' on the slot `bycolumn` set to `TRUE` or `FALSE`.
#' 
#' ## Methods
#' There are a couple of methods that intend to simplify the handling of data 
#' for the user. These methods are listed below. 
#' 
#' ### Show
#' * `show()` gives a short summary of the object's slots.
#' 
#' ### Getters
#' * `getY()` returns the `y` slot.
#' * `getColY()` returns the `y` slot in column format independent of 
#'   `bycolumn`.
#' * `getRowY()` returns the `y` slot in row format independent of `bycolumn`.
#' * `getN()` returns the `N` slot.
#' * `getr()` returns the `r` slot.
#' * `getS()` returns the `S` slot.
#' * `getColS()` returns the `S` slot in column format independent of 
#'   `bycolumn`.
#' * `getRowS()` returns the `S` slot in row format independent of `bycolumn`.
#' * `getBycolumn()` returns the `bycolumn` slot.
#' * `getName()` returns the `name` slot.
#' * `getType()` returns the `type` slot.
#' * `getSim()` returns the `sim` slot.
#' * `getExp()` returns the `exp` slot.
#' * `getColExp()` returns the `y` slot in column format independent of 
#'   `bycolumn`.
#' * `getRowExp()` returns the `y` slot in row format independent of `bycolumn`.
#' * `getT()` returns the `T` slot.
#' * `getColT()` returns the `T` slot in column format independent of 
#'   `bycolumn`.
#' * `getRowT()` returns the `T` slot in row format independent of `bycolumn`.
#' 
#' ### Setters
#' All setters help the user to set the slots in the right format and with the 
#' correct class (integer, matrix, etc.). It is internally checked, if the 
#' new value fits the other slots of the object. 
#' 
#' * `setY()<-` sets the `y` slot.
#' * `setN()<-` sets the `N` slot.
#' * `setR()<-` sets the `r` slot.
#' * `setS()<-` sets the `S` slot.
#' * `setBycolumn` sets the `bycolumn` slot.
#' * `setName()<-` sets the `name` slot.
#' * `setType()<-` sets the `type` slot.
#' * `setSim()<-` sets the `sim` slot.
#' * `setExp()<-` sets the `exp` slot.
#' * `setT()<-` sets the `T` slot.
#' 
#' ### Checking methods
#' The checking methods are provided to allow a user to integrate the `finmix` 
#' classes more easily into a larger code basis. They check, if the slots are 
#' available and return a `logical`. 
#'  
#' * `hasY()` checks, if slot `y` is not empty.
#' * `hasS()` checks, if slot `S` is not empty.
#' * `hasExp()` checks, if the slot `exp` is not empty.
#' * `hasT()` checks, if the slot `T` is not empty.
#' 
#' ### Plotting
#' The plotting function should help the user to get an impression of how the 
#' data in the `fdata` object is distributed. This is important for evaluating 
#' what kind of distribution to choose and how many mixture components to test 
#' for.
#' 
#' * `plot(x, y, dev=TRUE, ...)` plots the observations in the `y` slot. If the 
#'   `type` is `"discrete"` a [barplot()] is shown. In the `"continuous"` case 
#'   the plot depends on the number of dimensions: if the dimension `r` of the 
#'   data is one, a [histogram()] shows the distribution of the observations. 
#'   In case of a two-dimensional data set, histograms of the marginal 
#'   distributions are plotted together with a scatter [plot()] and a 
#'   two-dimensional kernel-density (see [bkde2D()]). In case of a multivariate 
#'   data set with more than two dimensions a [pairs()] plot is returned. The 
#'   argument `dev` should be put to `FALSE` if the output should be in a file. 
#'   `...` allows the user to pass further arguments to the internal functions.
#' 
#' @slot y A matrix containing the observations for finite mixture estimation. 
#'   Can be by column or row depending on the slot `bycolumn`.
#' @slot N An integer holding the number of observations.
#' @slot r An integer defining the dimension of the data. Only for multivariate
#'   distributions like `normult` or `studmult` the dimension is 
#'   larger one. 
#' @slot S A matrix containing the indicators of the data. If the `fdata` class
#'   contains indicators estimation is performed with a fixed indicator 
#'   approach.
#' @slot bycolumn A logical indicating if the data in `y` and `S` is sorted by
#'   by column (`TRUE`) or row (`FALSE`).
#' @slot name A character specifying a name for the data. Optional.
#' @slot type A character specifying the data type: either `discrete` for 
#'   discrete data or `continuous` for continuous data. The two data types are
#'   treated differently when calculating data moments. 
#' @slot sim A logical indicating, if the data was simulated. 
#' @slot exp A matrix containing the *exposures* of Poisson data.
#' @slot T A matrix containing the (optional) repetitions of binomial or Poisson
#'   data. Must be of type integer. 
#' @exportClass fdata
#' @rdname fdata-class
#' 
#' @seealso 
#' * [fdata()] for the class constructor
#' * [model-class] for the class from which data can be simulated
#' * [simulate()] for the method of the `model` class simulating data from a 
#'   finite mixture model
.fdata <- setClass("fdata",
  representation(
    y = "matrix",
    N = "integer",
    r = "integer",
    S = "matrix",
    bycolumn = "logical",
    name = "character",
    type = "character",
    sim = "logical",
    exp = "matrix",
    T = "matrix"
  ),
  validity = function(object) {
    .valid.Fdata(object)
    ## else: ok
    TRUE
  },
  prototype(
    y = matrix(),
    N = integer(),
    r = integer(),
    S = matrix(),
    bycolumn = logical(),
    name = character(),
    type = "discrete",
    sim = logical(),
    exp = matrix(),
    T = matrix()
  )
)

## Constructor for the data class ##
#' Constructs an `fdata` object
#' 
#' @description 
#' Calling [fdata()] constructs an `fdata` object. Can be called without 
#' arguments.
#' 
#' @param y A matrix containing the observations for finite mixture estimation. 
#'   Can be by column or row depending on the slot `bycolumn`. 
#' @param N An integer holding the number of observations.
#' @param r An integer defining the dimension of the data. Only for multivariate
#'   distributions like `normult` or `studmult` the dimension is 
#'   larger one. 
#' @param S A matrix containing the indicators of the data. If the `fdata` class
#'   contains indicators estimation is performed with a fixed indicator 
#'   approach.
#' @param bycolumn A logical indicating if the data in `y` and `S` is sorted by
#'   by column (`TRUE`) or row (`FALSE`).
#' @param name A character specifying a name for the data. Optional.
#' @param type A character specifying the data type: either `discrete` for 
#'   discrete data or `continuous` for continuous data. The two data types are
#'   treated differently when calculating data moments. 
#' @param sim A logical indicating, if the data was simulated. 
#' @param exp A matrix containing the *exposures* of Poisson data.
#' @param T A matrix containing the (optional) repetitions of binomial or Poisson
#'   data. Must be of type integer. 
#' @export
#' 
#' @examples
#' # Call the constructor without arguments.
#' f_data <- fdata()
#' 
#' # Create simulated data.
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' 
#' @seealso [fdata] class that describes the slots and the getters, setters and
#'   and checkers
"fdata" <- function(y = matrix(), N = 1, r = 1, S = matrix(),
                    bycolumn = TRUE, name = character(),
                    type = "discrete", sim = FALSE,
                    exp = matrix(), T = matrix()) {
  y <- as.matrix(y)
  .check.y.Fdata(y)
  if (missing(type)) {
    type <- .check.type.Fdata(y)
  }
  if (missing(bycolumn)) {
    bycolumn <- .check.bycolumn.Fdata(y, S, exp, T)
  }
  if (missing(N)) {
    N <- .check.N.Fdata(y, S, exp, T, bycolumn)
  } else {
    N <- as.integer(N)
  }
  if (missing(r)) {
    r <- .check.r.Fdata(y, bycolumn)
  } else {
    r <- as.integer(r)
  }
  if (!all(is.na(S))) {
    S <- .check.S.Fdata(S)
  } else {
    storage.mode(S) <- "integer"
  }
  if (!all(is.na(exp))) {
    exp <- .check.exp.Fdata(exp)
  }
  if (!all(is.na(T))) {
    T <- .check.T.Fdata(T)
  } else {
    storage.mode(T) <- "integer"
  }
  .fdata(
    y = y, N = N, r = r, S = S,
    bycolumn = bycolumn, name = name,
    type = type, sim = sim, exp = exp,
    T = T
  )
}

#' Plots the data
#' 
#' @description
#' [plot()] plots the data in an [fdata] object by either a histogram in case of
#' continuous data or a barplot in case of discrete data.
#' 
#' @param x An `fdata` object. Cannot be empty.
#' @param y Unused.
#' @param dev A logical indicating if the plot should be output via a graphical
#'   device.
#' @param ... Further arguments passed to the plotting functions `hist` or 
#'   `barplot`. 
#' @exportMethod plot
#' @keywords internal
#' 
#' @examples
#' # Generate Poisson data and plot it. 
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' plot(f_data)
#' 
#' @seealso [fdata] class
setMethod(
  "plot", signature(
    x = "fdata",
    y = "missing"
  ),
  function(x, y, dev = TRUE, ...) {
    hasY(x, verbose = TRUE)
    if (.check.grDevice() && dev) {
      dev.new(title = "Histograms")
    }
    if (x@type == "discrete") {
      .plot.discrete.Fdata(x)
    } else { ## continuous
      .plot.continuous.Fdata(x, dev)
    }
  }
)

#' Shows a summary of an `fdata` object.
#' 
#' Calling [show()] on an `fdata` object gives an overview of the different 
#' slots and dimensions. 
#' 
#' @param object An `fdata` object.
#' @returns A console output listing the slots and summary information about
#'   each of them. 
#' @exportMethod show
#' @keywords internal
#' 
#' @examples 
#' # Generate some Poisson data and show the `fdata` object
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' show(f_data)
#' 
#' @seealso [fdata] class for an overview of the slots
setMethod(
  "show", "fdata",
  function(object) {
    name <- ifelse(length(object@name) == 0, "fdata",
      object@name
    )
    cat("Object '", name, "'\n", sep = "")
    cat("     class       :", class(object), "\n")
    cat(
      "     y           :",
      paste(dim(object@y), collapse = "x"), "\n"
    )
    cat("     bycolumn    :", object@bycolumn, "\n")
    cat("     N           :", object@N, "\n")
    cat("     r           :", object@r, "\n")
    if (hasS(object)) {
      cat(
        "     S           :",
        paste(dim(object@S), collapse = "x"), "\n"
      )
    }
    cat("     type        :", object@type, "\n")
    cat("     sim         :", object@sim, "\n")
    if (hasExp(object)) {
      cat(
        "     exp         :",
        paste(dim(object@exp), collapse = "x"), "\n"
      )
    }
    if (hasT(object)) {
      cat(
        "     T           :",
        paste(dim(object@T), collapse = "x"), "\n"
      )
    }
  }
)

#' Checker method for `y` slot of an `fdata` object. 
#' 
#' @description 
#' `hasY()` checks, if the object contains `y` data.
#' 
#' @param object An `fdata` object. 
#' @param verbose A logical indicating, if the function should print out 
#'   messages.
#' @returns Either `FALSE`/`TRUE`, if `verbose` is `FALSE` and the `y` slot is 
#'   empty or filled or a message, if `verbose` is `TRUE`.
#' @exportMethod hasY
#' @keywords internal
#' 
#' @examples 
#' # Generate an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' hasY(f_data)
#' 
#' @seealso [fdata] class for an overview of its slots
setMethod(
  "hasY", "fdata",
  function(object, verbose = FALSE) {
    if (!all(is.na(object@y))) {
      return(TRUE)
    } else {
      if (verbose) {
        stop(paste("Slot 'y' in 'fdata' object ",
          "is empty.",
          sep = ""
        ))
      } else {
        return(FALSE)
      }
    }
  }
)

#' Checker method for `S` slot of an `fdata` object. 
#' 
#' @description 
#' `hasS()` checks, if the object contains `S` data.
#' 
#' @param object An `fdata` object. 
#' @param verbose A logical indicating, if the function should print out 
#'   messages.
#' @returns Either `FALSE`/`TRUE`, if `verbose` is `FALSE` and the `S` slot is 
#'   empty or filled or a message, if `verbose` is `TRUE`.
#' @exportMethod hasS
#' @keywords internal 
#' 
#' @examples 
#' # Generate an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' hasS(f_data)
#' 
#' @seealso 
#' * [fdata-class] for the class definition
setMethod(
  "hasS", "fdata",
  function(object, verbose = FALSE) {
    if (!all(is.na(object@S))) {
      return(TRUE)
    } else {
      if (verbose) {
        stop(paste("Slot 'S' in 'fdata' object ",
          "is empty.",
          sep = ""
        ))
      } else {
        return(FALSE)
      }
    }
  }
)

#' Checker method for `exp` slot of an `fdata` object. 
#' 
#' @description 
#' `hasExp()` checks, if the object contains `exp` data.
#' 
#' @param object An `fdata` object. 
#' @param verbose A logical indicating, if the function should print out 
#'   messages.
#' @returns Either `FALSE`/`TRUE`, if `verbose` is `FALSE` and the `exp` slot is 
#'   empty or filled or a message, if `verbose` is `TRUE`.
#' @exportMethod hasExp
#' @keywords internal
#' 
#' @examples 
#' # Generate an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' hasExp(f_data)
#' 
#' @seealso 
#' * [fdata-class] for the class definition
setMethod(
  "hasExp", "fdata",
  function(object, verbose = FALSE) {
    if (!all(is.na(object@exp))) {
      return(TRUE)
    } else {
      if (verbose) {
        stop(paste("Slot 'exp' in 'fdata' object ",
          "is empty.",
          sep = ""
        ))
      } else {
        return(FALSE)
      }
    }
  }
)

#' Checker method for `T` slot of an `fdata` object. 
#' 
#' @description 
#' `hasY()` checks, if the object contains `T` data.
#' 
#' @param object An `fdata` object. 
#' @param verbose A logical indicating, if the function should print out 
#'   messages.
#' @returns Either `FALSE`/`TRUE`, if `verbose` is `FALSE` and the `T` slot is 
#'   empty or filled or a message, if `verbose` is `TRUE`.
#' @exportMethod hasT
#' @keywords internal
#' @examples 
#' # Generate an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' hasT(f_data)
#' 
#' @seealso 
#' * [fdata-class] for the class defintion
setMethod(
  "hasT", "fdata",
  function(object, verbose = FALSE) {
    if (!all(is.na(object@T))) {
      return(TRUE)
    } else {
      if (verbose) {
        stop(paste("Slot @T in 'fdata' object ",
          "is empty.",
          sep = ""
        ))
      } else {
        return(FALSE)
      }
    }
  }
)

### getCol/getRow: These methods return the data in the slots @y,
### @S, @exp and @T either as column-ordered or row-ordered matrix.
#' Getter method of `fdata` class.
#' 
#' Returns the `y` slot as a column-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `y` slot of the `object` as a column-ordered matrix.
#' @exportMethod getColY
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getColY(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getColY", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(object@y)
    } else {
      return(t(object@y))
    }
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `y` slot as a row-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `y` slot of the `object` as a row-ordered matrix.
#' @exportMethod getRowY
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getRowY(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getRowY", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(t(object@y))
    } else {
      return(object@y)
    }
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `S` slot as a column-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `S` slot of the `object` as a column-ordered matrix.
#' @exportMethod getColS
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getColS(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getColS", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(object@S)
    } else {
      return(t(object@S))
    }
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `S` slot as a row-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `S` slot of the `object` as a row-ordered matrix.
#' @exportMethod getRowS
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getRowS(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getRowS", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(t(object@S))
    } else {
      return(object@S)
    }
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `exp` slot as a column-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `exp` slot of the `object` as a column-ordered matrix.
#' @exportMethod getColExp
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getColExp(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getColExp", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(object@exp)
    } else {
      return(t(object@exp))
    }
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `exp` slot as a row-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `exp` slot of the `object` as a row-ordered matrix.
#' @exportMethod getRowExp
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getRowExp(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getRowExp", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(t(object@exp))
    } else {
      return(object@exp)
    }
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `T` slot as a column-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `T` slot of the `object` as a column-ordered matrix.
#' @exportMethod getColT
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getColT(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getColT", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(object@T)
    } else {
      return(t(object@T))
    }
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `T` slot as a row-ordered matrix. 
#' 
#' @param object An `fdata` object.
#' @returns The `T` slot of the `object` as a row-ordered matrix.
#' @exportMethod getRowT
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getRowT(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getRowT", "fdata",
  function(object) {
    if (object@bycolumn) {
      return(t(object@T))
    } else {
      return(object@T)
    }
  }
)

## Setters and Getters as a user interface to manipulate the slots
## Combined Getter and Setter
#' Getter method of `fdata` class.
#' 
#' Returns the `y` slot in the order defined by the slot `bycolumn`. 
#' 
#' @param object An `fdata` object.
#' @returns The `y` slot of the `object` in the order defined `bycolumn`.
#' @exportMethod getY
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getY(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getY", "fdata",
  function(object) {
    return(object@y)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `N` slot of an `fdata` object. 
#' 
#' @param object An `fdata` object.
#' @returns The `N` slot of the `object`.
#' @exportMethod getN
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getN(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getN", "fdata",
  function(object) {
    return(object@N)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `r` slot of an `fdata` object. 
#' 
#' @param object An `fdata` object.
#' @returns The `r` slot of the `object`.
#' @exportMethod getR
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getR(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getR", "fdata",
  function(object) {
    return(object@r)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `S` slot in the order defined by the slot `bycolumn`. 
#' 
#' @param object An `fdata` object.
#' @returns The `S` slot of the `object` in the order defined `bycolumn`.
#' @exportMethod getS
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getS(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getS", "fdata",
  function(object) {
    return(object@S)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `bycolumn` slot of an `fdata` object. 
#' 
#' @param object An `fdata` object.
#' @returns The `bycolumn` slot of the `object`.
#' @exportMethod getBycolumn
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getBycolumn(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getBycolumn", "fdata",
  function(object) {
    return(object@bycolumn)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `name` slot of an `fdata` object. 
#' 
#' @param object An `fdata` object.
#' @returns The `name` slot of the `object`.
#' @exportMethod getName
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getName(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getName", "fdata",
  function(object) {
    return(object@name)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `type` slot of an `fdata` object. 
#' 
#' @param object An `fdata` object.
#' @returns The `type` slot of the `object`.
#' @exportMethod getType
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getType(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getType", "fdata",
  function(object) {
    return(object@type)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `sim` slot of an `fdata` object. 
#' 
#' @param object An `fdata` object.
#' @returns The `sim` slot of the `object`.
#' @exportMethod getSim
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getSim(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getSim", "fdata",
  function(object) {
    return(object@sim)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `exp` slot in the order defined by the slot `bycolumn`. 
#' 
#' @param object An `fdata` object.
#' @returns The `exp` slot of the `object` in the order defined `bycolumn`.
#' @exportMethod getExp
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getExp(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getExp", "fdata",
  function(object) {
    return(object@exp)
  }
)

#' Getter method of `fdata` class.
#' 
#' Returns the `T` slot in the order defined by the slot `bycolumn`. 
#' 
#' @param object An `fdata` object.
#' @returns The `T` slot of the `object` in the order defined `bycolumn`.
#' @exportMethod getT
#' @keywords internal
#' 
#' @examples 
#' # Create an fdata object with Poisson data
#' f_data <- fdata(y = rpois(100, 312), sim = TRUE)
#' getT(f_data)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setMethod(
  "getT", "fdata",
  function(object) {
    return(object@T)
  }
)

## Setters ##
#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `y` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `y` should be set.
#' @param value A matrix that should be set as `y` slot of the `fdata` object.
#' @returns The `fdata` object with slot `y` set to `value` or an error message
#'   if the `value` cannot be set as slot `y`.
#' @exportMethod setY<-
#' @keywords internal
#' 
#' @examples
#' f_data <- fdata()
#' setY(f_data) <- rpois(100, 312)
#' 
#' @seealso [fdata-class] for all slots of the `fdata` class
setReplaceMethod(
  "setY", "fdata",
  function(object, value) {
    value <- as.matrix(value)
    .check.y.Fdata(value)
    if (object@bycolumn && NROW(value) < NCOL(value)) {
      object@y <- t(value)
    } else {
      object@y <- value
    }
    if (object@bycolumn) {
      object@N <- NROW(object@y)
      object@r <- NCOL(object@y)
    } else {
      object@N <- NCOL(object@y)
      object@r <- NROW(object@y)
    }
    .init.valid.Fdata(object)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `N` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `N` should be set.
#' @param value An integer that should be set as `N` slot of the `fdata` object.
#' @returns The `fdata` object with slot `N` set to `value` or an error message
#'   if the `value` cannot be set as slot `N`.
#' @exportMethod setN<-
#' @keywords internal
#' 
#' @examples
#' f_data <- fdata()
#' setN(f_data) <- as.integer(100)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setN", "fdata",
  function(object, value) {
    object@N <- as.integer(value)
    .init.valid.Fdata(object)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `R` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `R` should be set.
#' @param value An integer that should be set as `R` slot of the `fdata` object.
#' @returns The `fdata` object with slot `R` set to `value` or an error message
#'   if the `value` cannot be set as slot `R`.
#' @exportMethod setR<-
#' @keywords internal
#' 
#' @examples
#' f_data <- fdata()
#' setR(f_data) <- as.integer(2)
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setR", "fdata",
  function(object, value) {
    object@r <- as.integer(value)
    .init.valid.Fdata(object)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `S` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` object, whose slot `S` should be set.
#' @param value A matrix that should be set as `S` slot of the `fdata` object. 
#'   Has to be of type integer.
#' @returns The `fdata` object with slot `S` set to `value` or an error message
#'   if the `value` cannot be set as slot `S`.
#' @exportMethod setS<-
#' @keywords internal
#' 
#' @examples
#' # Generate an empty fdata object.
#' f_data <- fdata()
#' setS(f_data) <- matrix(sample.int(4, 100, replace = TRUE))
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setS", "fdata",
  function(object, value) {
    value <- .check.S.Fdata(value)
    if (object@bycolumn && NROW(value) > NCOL(value)) {
      object@S <- value
    } else {
      object@S <- t(value)
    }
    .init.valid.Fdata(object)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `bycolumn` of an `fdata` object and validates the slot data 
#' before setting. 
#' 
#' @param object An `fdata` objects, whose slot `bycolumn` should be set.
#' @param value A logical that should be set as `bycolumn` slot of the `fdata` 
#'   object. 
#' @returns The `fdata` object with slot `bycolumn` set to `value` or an error message
#'   if the `value` cannot be set as slot `bycolumn`.
#' @exportMethod setBycolumn<-
#' @keywords internal
#' 
#' @examples
#' # Generate an empty fdata object.
#' f_data <- fdata()
#' setBycolumn(f_data) <- TRUE
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setBycolumn", "fdata",
  function(object, value) {
    .check.setBycolumn.Fdata(value)
    tmp.bycolumn <- object@bycolumn
    if (tmp.bycolumn != value) {
      object@bycolumn <- value
      if (!all(is.na(object@y))) {
        object@y <- t(object@y)
      }
      if (!all(is.na(object@S))) {
        object@S <- t(object@S)
      }
      if (!all(is.na(object@exp))) {
        object@exp <- t(object@exp)
      }
      if (!all(is.na(object@T))) {
        object@T <- t(object@T)
      }
    }
    .init.valid.Fdata(object)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `name` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `name` should be set.
#' @param value A matrix that should be set as `name` slot of the `fdata` object. 
#'   Has to be of type integer.
#' @returns The `fdata` object with slot `name` set to `value` or an error message
#'   if the `value` cannot be set as slot `name`.
#' @exportMethod setName<-
#' @keywords internal
#' 
#' @examples
#' # Generate an empty fdata object.
#' f_data <- fdata()
#' setName(f_data) <- "poisson_data"
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setName", "fdata",
  function(object, value) {
    object@name <- as.character(value)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `type` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `type` should be set.
#' @param value A character that should be set as `type` slot of the `fdata` object. 
#'   Has to be of type integer.
#' @returns The `fdata` object with slot `type` set to `value` or an error message
#'   if the `value` cannot be set as slot `type`.
#' @exportMethod setType<-
#' @keywords internal
#' 
#' @examples
#' # Generate an empty fdata object.
#' f_data <- fdata()
#' setType(f_data) <- "discrete"
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setType", "fdata",
  function(object, value) {
    object@type <- as.character(value)
    .valid.type.Fdata(object)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `sim` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `sim` should be set.
#' @param value A logical that should be set as `sim` slot of the `fdata` object. 
#'   Has to be of type integer.
#' @returns The `fdata` object with slot `sim` set to `value` or an error message
#'   if the `value` cannot be set as slot `sim`.
#' @exportMethod setSim<-
#' @keywords internal
#' 
#' @examples
#' # Generate an empty fdata object.
#' f_data <- fdata()
#' setSim(f_data) <- TRUE
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setSim", "fdata",
  function(object, value) {
    .check.setSim.Fdata(value)
    object@sim <- value
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `exp` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `exp` should be set.
#' @param value A matrix that should be set as `exp` slot of the `fdata` object. 
#'   Has to be of type integer.
#' @returns The `fdata` object with slot `exp` set to `value` or an error message
#'   if the `value` cannot be set as slot `exp`.
#' @exportMethod setExp<-
#' @keywords internal
#' 
#' @examples
#' # Generate an empty fdata object.
#' f_data <- fdata()
#' setExp(f_data) <- matrix(rep(100, 100))
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setExp", "fdata",
  function(object, value) {
    value <- matrix(value)
    value <- .check.exp.Fdata(value)
    if (object@bycolumn && NROW(value) > NCOL(value)) {
      object@exp <- value
    } else {
      object@exp <- t(value)
    }
    .init.valid.Fdata(object)
    return(object)
  }
)

#' Setter method of `fdata` class
#' 
#' @description
#' Sets the slot `T` of an `fdata` object and validates the slot data before 
#' setting. 
#' 
#' @param object An `fdata` objects, whose slot `T` should be set.
#' @param value A matrix that should be set as `T` slot of the `fdata` object. 
#'   Has to be of type integer.
#' @returns The `fdata` object with slot `T` set to `value` or an error message
#'   if the `value` cannot be set as slot `T`.
#' @exportMethod setT<-
#' @keywords internal
#' 
#' @examples
#' # Generate an empty fdata object.
#' f_data <- fdata()
#' setT(f_data) <- matrix(rep(100, 100))
#' 
#' @seealso [fdata] for all slots of the `fdata` class
setReplaceMethod(
  "setT", "fdata",
  function(object, value) {
    value <- matrix(value)
    value <- .check.T.Fdata(value)
    if (object@bycolumn && NROW(value) > NCOL(value)) {
      object@T <- value
    } else {
      object@T <- t(value)
    }
    .init.valid.Fdata(object)
    return(object)
  }
)

### Private functions.
### These functions are not exported.

### Checking.
### Check data: The data @y has to either of type 'integer'
### or of type 'numeric'.
#' Checks validity of slot `y`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.
#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setY()] the setter of the slot `y`.
".check.y.Fdata" <- function(y) {
  if (!all(is.na(y))) {
    ## Only data of type 'numeric' or
    ## 'integer' is accepted.
    if (!is.numeric(y) && !is.integer(y)) {
      stop(paste("Argument 'y' must be of type ",
        "'numeric' or 'integer'.",
        sep = ""
      ))
    }
  }
}

### Check type: The type @type has to be either 'discrete' or
### 'continuous'. If @y is of storage mode 'integer' @type
### is set to 'discrete', else 'continuous'.
### If @y is NA for all entries, the type is the default:
### 'discrete'.
#' Checks validity of slot `type`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.
#'   
#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setType()] the setter of the slot `type`.
".check.type.Fdata" <- function(y) {
  if (!all(is.na(y))) {
    if (is.integer(y)) {
      return("discrete")
    } else {
      return("continuous")
    }
  } else {
    return("discrete")
  }
}

### Check bycolumn: The data is stored either by row or by column.
### If the data in @y has more rows than columns, it is assumed,
### that it must be stored by column. Otherwise, it is assumed,
### that it must be stored by row.
### If rows are equal or less than columns, @bycolumn is set to
### FALSE.
### If the data in @y is empty, it is checked in the same way if
### bycolumn can be derived from @S, @exp or @T. If any data slot is
### emtpy the default is used: TRUE.
#' Checks validity of slot `bycolumn`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.

#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setBycolumn()] the setter of the slot `bycolumn`.
".check.bycolumn.Fdata" <- function(y, S, exp, T) {
  if (!all(is.na(y))) {
    if (NROW(y) > NCOL(y)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    if (!all(is.na(S))) {
      if (NROW(S) > NCOL(S)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      if (!all(is.na(exp))) {
        if (NROW(exp) > NCOL(exp)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      } else {
        if (!all(is.na(T))) {
          if (NROW(T) > NCOL(T)) {
            return(TRUE)
          } else {
            return(FALSE)
          }
        } else {
          return(TRUE)
        }
      }
    }
  }
}

### Check N: The number of observations @N of the dataset @y is
### set after @bycolumn. So, if @bycolumn is TRUE, the rows are
### assumed to be the number of observations. Otherwise, columns
### of @y are assumed to define @N.
#' Checks validity of slot `N`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.

#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setN()] the setter of the slot `N`.
".check.N.Fdata" <- function(y, S, exp, T, bycolumn) {
  if (!all(is.na(y))) {
    if (bycolumn) {
      NROW(y)
    } else {
      NCOL(y)
    }
  } else {
    if (!all(is.na(S))) {
      if (bycolumn) {
        NROW(S)
      } else {
        NCOL(S)
      }
    } else {
      if (!all(is.na(exp))) {
        if (bycolumn) {
          NROW(exp)
        } else {
          NCOL(exp)
        }
      } else {
        if (!all(is.na(T))) {
          if (bycolumn) {
            NROW(T)
          } else {
            NCOL(T)
          }
        } else {
          return(as.integer(1))
        }
      }
    }
  }
}

### Check r: The number @r of variables in the dataset @y is set
### after @bycolumn. So, if @bycolum is TRUE, the columns are assumed
### to determine the number of variables. Otherwise, rows of @y are
### assumed to define @r.
#' Checks validity of slot `r`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.

#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setR()] the setter of the slot `r`.
".check.r.Fdata" <- function(y, bycolumn) {
  if (!all(is.na(y))) {
    if (bycolumn) {
      NCOL(y)
    } else {
      NROW(y)
    }
  } else {
    return(as.integer(1))
  }
}

### Check S: Indicators must be of type 'integer'. If this is the case
### the indicators are turned into a matrix object with storage mode
### 'integer'.
#' Checks validity of slot `S`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.

#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setS()] the setter of the slot `S`.
".check.S.Fdata" <- function(S) {
  if (!all(is.na(S))) {
    if (!is.numeric(S)) {
      stop(paste("Wrong type of slot 'S' in 'fdata' object. ",
        "Indicators must be of type 'integer'.",
        sep = ""
      ))
    } else {
      S <- as.matrix(S)
      storage.mode(S) <- "integer"
      return(S)
    }
  } else {
    return(S)
  }
}

### Check T: Repetitions must be of type 'integer'. If this is the case
### the repetitions are turned into a matrix object with storage mode
### 'integer'.
#' Checks validity of slot `T`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.

#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setT()] the setter of the slot `T`.
".check.T.Fdata" <- function(T) {
  if (!all(is.na(T))) {
    if (!is.numeric(T)) {
      stop(paste("Wrong type of slot 'T' in 'fdata' object. ",
        "Repetitions must be of type 'integer'.",
        sep = ""
      ))
    } else {
      T <- as.matrix(T)
      storage.mode(T) <- "integer"
      return(T)
    }
  }
}

### Check exp: Exposures must be of of type 'numeric'. If this is
### the case exposures are turned into a matrix.
#' Checks validity of slot `exp`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.
#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setY()] the setter of the slot `exp`.
".check.exp.Fdata" <- function(exp) {
  if (!all(is.na(exp))) {
    if (!is.numeric(exp)) {
      stop(paste("Wrong type of slot 'exp' in 'fdata' object. ",
        "Exposures must be of type 'numeric'.",
        sep = ""
      ))
    } else {
      exp <- as.matrix(exp)
      return(exp)
    }
  }
}

### Check bycolumn: @bycolumn has to be of type 'logical'. If this is not
### the case an error is thrown.
#' Checks validity of slot `bycolumn`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.
#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setBycolumn()] the setter of the slot `bycolumn`.
".check.setBycolumn.Fdata" <- function(value) {
  if (!is.logical(value)) {
    stop(paste("Wrong specification of value for slot 'bycolumn' ",
      "in 'fdata' object. 'bycolumn' must be of type ",
      "'logical'.",
      sep = ""
    ))
  }
}

### Check sim: @sim has to be of type 'logical'. If this is not
### the case an error is thrown.
#' Checks validity of slot `sim`
#' 
#' @param y An object passed in by the user. 
#' @returns None. Checks for validity and if validity is not ensured throws an
#'   error.
#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [setSim()] the setter of the slot `sim`.
".check.setSim.Fdata" <- function(value) {
  if (!is.logical(value)) {
    stop(paste("Wrong specification of value for slot 'sim' ",
      "in 'fdata' object. 'sim' must be of type ",
      "'logical'.",
      sep = ""
    ))
  }
}

### Plot functions
### Discrete data: Only univariate discrete data is
### is implemented. The functions plots a barplot.
### If the data in @y has names given, these names
### are used in the plot.
#' Plots discrete data of an `fdata` object
#' 
#' @param obj An `fdata` object. Must contain data.
#' @returns A barplot. 
#' 
#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [plot()] the plot function of the `fdata` class.
#' * [barplot()] the default plotting function for bar plots in R. 
".plot.discrete.Fdata" <- function(obj) {
  if (hasY(obj, verbose = TRUE)) {
    datam <- getColY(obj)
  }
  if (hasExp(obj)) {
    exp <- getColExp(obj)
    datam <- datam * exp
  }
  barplot(table(datam),
    col = "gray65",
    border = "white", cex = 0.7,
    cex.axis = 0.7, xlab = "", main = "",
    cex.lab = 0.7
  )
  if (!is.null(colnames(datam))) {
    col.names <- colnames(datam)
  } else {
    col.names <- c("")
  }
  mtext(
    side = 1, col.names, cex = 0.7,
    line = 3
  )
}

### Continuous data: Either the data is one-dimensional or
### multi-dimensional.
### In the one-dimensional case a histogram of the data is
### plotted.
### In the two-dimensional case a bivariate kernel density
### estimation is used to return a contour plot and a
### perspective plot of the density.
### In the case of higher-dimensional data, the functions
### returns histograms for all variables in @y and a pairs
### diagram: a matrix containing scatter plots for all
### variables' combinations.
#' Plots discrete data of an `fdata` object
#' 
#' @description 
#' Continuous data: Either the data is one-dimensional or multi-dimensional. In 
#' the one-dimensional case a histogram of the data is plotted. In the 
#' two-dimensional case a bivariate kernel density estimation is used to return 
#' a contour plot and a perspective plot of the density. In the case of 
#' higher-dimensional data, the functions returns histograms for all variables 
#' in `@@y` and a pairs diagram: a matrix containing scatter plots for all 
#' variables' combinations.
#' 
#' @param obj An `fdata` object. Must contain data.
#' @returns A histogram. 
#' @importFrom KernSmooth bkde2D
#' @importFrom stats sd
#' @noRd
#' 
#' @seealso 
#' * [fdata()] the constructor of the `fdata` class.
#' * [plot()] the plot function of the `fdata` class.
#' * [hist()] the default plotting function for histogram plots in R. 
".plot.continuous.Fdata" <- function(obj, dev) {
  datam <- getColY(obj)
  if (is.null(colnames(datam))) {
    colnames(datam) <- paste("y", seq(1, ncol(datam)),sep = "")
  }
  if (obj@r == 1) {
    .symmetric.Hist(datam, colnames(datam))
  } else if (obj@r == 2) { ## 2-dimensional
    .symmetric.Hist(datam, colnames(datam))
    if (.check.grDevice() && dev) {
      dev.new(title = "Contour plot")
    }
    par(
      mfrow = c(1, 2), mar = c(2, 2, 2, 3),
      oma = c(4, 5, 1, 5)
    )
    plot(datam[, 1], datam[, 2],
      col = "gray47",
      cex = 0.7, cex.axis = 0.7,
      pch = 20, xlab = "", ylab = "",
      main = ""
    )
    mtext(
      side = 1, colnames(datam)[1],
      cex = 0.7, line = 3
    )
    mtext(
      side = 2, colnames(datam)[2],
      cex = 0.7, line = 3
    )
    d <- bkde2D(datam,
      bandwidth = c(
        sd(datam[, 1]),
        sd(datam[, 2])
      )
    )
    contour(d$x1, d$x2, d$fhat,
      col = "gray47",
      cex = 0.7, cex.axis = 0.7,
      xlab = "", ylab = ""
    )
    mtext(
      side = 1, colnames(datam)[1],
      cex = 0.7, line = 3
    )
    mtext(
      side = 2, colnames(datam)[2],
      cex = 0.7, line = 3
    )
    if (.check.grDevice() && dev) {
      dev.new("Perspective plot")
      par(mfrow = c(1, 1))
    }
    if (!is.null(colnames(datam))) {
      col.names <- colnames(datam)
    } else {
      col.names <- c("", "")
    }
    persp(d$x1, d$x2, d$fhat,
      main = "",
      xlab = col.names[1], ylab = col.names[2],
      zlab = "", col = "gray65",
      border = "gray47", theta = 55, phi = 30,
      expand = 1.6, lphi = 190, ltheta = 90,
      r = 40, d = 0.1, cex = 0.7, cex.axis = 0.7,
      cex.lab = 0.7, ticktype = "detailed"
    )
  } else { ## multivariate distribution
    .symmetric.Hist(datam, colnames(datam))
    if (.check.grDevice() && dev) {
      dev.new(title = "Pairs")
    }
    pairs(datam,
      col = "gray47", pch = 20,
      cex = 0.7, cex.axis = 0.7, cex.labels = 1.3
    )
  }
}

### Validity
### Initial data: If the 'fdata' object is modified via setters,
### the user may define the slots step by step.
### 'fdata()'. To avoid cumbersome behavior of slot setting,
### only warnings are thrown.
#' Checks consistency of `fdata` object
#' 
#' @description
#' Checks the consistency of `fdata` object when the constructor or a setter is
#' called. The different slots are strongly related and it has to be checked if 
#' the setting of one slot does not interfere with the definition of another 
#' one. Is called during initialization of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [fdata()] for the constructor and the conditions for slots. 
".init.valid.Fdata" <- function(obj) {
  .init.valid.y.Fdata(obj)
  .init.valid.S.Fdata(obj)
  .init.valid.exp.Fdata(obj)
  .init.valid.T.Fdata(obj)
  .valid.type.Fdata(obj)
}

### Valid data: If later data objects are used in functions, the functions
### usually need fully specified and consistent slots of an 'fdata'
### object. For this case, the 'fdata' object can be checked with
### errors thrown in case of inconsistency of slots.
### Furthermore, the validity check during initialisation relies on fully
### specified 'fdata' objects and checks consistency of slots strongly.
#' Checks consistency of `fdata` object
#' 
#' @description
#' Checks the consistency of `fdata` object when the constructor or a setter is
#' called. The different slots are strongly related and it has to be checked if 
#' the setting of one slot does not interfere with the definition of another 
#' one. Is called during by setters of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [fdata()] for the constructor and the conditions for slots. 
".valid.Fdata" <- function(obj) {
  .valid.y.Fdata(obj)
  .valid.S.Fdata(obj)
  .valid.exp.Fdata(obj)
  .valid.T.Fdata(obj)
  .valid.type.Fdata(obj)
}

### Valid y: Data in @y must be of type 'integer' or 'numeric'. Further,
### the number of observations @N, the dimension of observations @r
### and the ordering @bycolumn must be group-consistent.
#' Checks consistency of slot `y` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `y` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called during initialization of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [fdata()] for the constructor and the conditions for slots. 
".init.valid.y.Fdata" <- function(obj) {
  if (!all(is.na(obj@y))) {
    if (!is.numeric(obj@y) && !is.integer(obj@y)) {
      stop(paste("Wrong type of slot 'y' in 'fdata' object. ",
        "Data must be of type 'numeric' or 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (obj@N != nrow(obj@y)) {
        warning(paste("Slot 'N' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is TRUE",
          sep = ""
        ))
      }
      if (obj@r != ncol(obj@y)) {
        warning(paste("Slot 'r' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is TRUE",
          sep = ""
        ))
      }
    } else {
      if (obj@N != ncol(obj@y)) {
        warning(paste("Slot 'N' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is FALSE",
          sep = ""
        ))
      }
      if (obj@r != nrow(obj@y)) {
        warning(paste("Slot 'r' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is FALSE",
          sep = ""
        ))
      }
    }
  }
}

#' Checks consistency of slot `y` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `y` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called from the setter of slot `y` of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [setY()] for the constructor and the conditions for slots. 
".valid.y.Fdata" <- function(obj) {
  if (!all(is.na(obj@y))) {
    if (!is.numeric(obj@y) && !is.integer(obj@y)) {
      stop(paste("Wrong type of slot 'y' in 'fdata' object. ",
        "Data must be of type 'numeric' or 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (obj@N != nrow(obj@y)) {
        stop(paste("Slot 'N' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is TRUE",
          sep = ""
        ))
      }
      if (obj@r != ncol(obj@y)) {
        stop(paste("Slot 'r' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is TRUE",
          sep = ""
        ))
      }
    } else {
      if (obj@N != ncol(obj@y)) {
        stop(paste("Slot 'N' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is FALSE",
          sep = ""
        ))
      }
      if (obj@r != nrow(obj@y)) {
        stop(paste("Slot 'r' does not match the ",
          "dimension of slot 'y' in 'fdata' ",
          "object if slot 'bycolumn' is FALSE",
          sep = ""
        ))
      }
    }
  }
}

### Valid S: Indicators in @S must be of type 'integer'. Further,
### the number of observations in @y and @S must be consistent.
### Consistency must also hold in regard to ordering defined by
### @bycolumn.
### Indicators must be positive integers. If any element of @S
### is smaller than one, an error is thrown.
#' Checks consistency of slot `S` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `S` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called during initialization of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [fdata()] for the constructor and the conditions for slots. 
".init.valid.S.Fdata" <- function(obj) {
  if (!all(is.na(obj@S))) {
    if (!is.integer(obj@S)) {
      stop(paste("Wrong type of slot 'S' in 'fdata' object. ",
        "Indicators must be of type 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (!all(is.na(obj@y))) {
        if (nrow(obj@S) != nrow(obj@y)) {
          warning(paste("Dimension of slot 'S' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is TRUE.",
            sep = ""
          ))
        }
      }
      if (ncol(obj@S) > 1) {
        warning(paste("Wrong dimension of slot 'S' if ",
          "slot 'bycolumn' is TRUE. Indicators ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    } else {
      if (!all(is.na(obj@y))) {
        if (ncol(obj@S) != ncol(obj@y)) {
          warning(paste("Dimension of slot 'S' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is FALSE.",
            sep = ""
          ))
        }
      }
      if (nrow(obj@S) > 1) {
        warning(paste("Wrong dimension of slot 'S' if ",
          "slot 'bycolumn' is FALSE. Indicators ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    }
    if (any(obj@S < 1)) {
      stop(paste("Wrong speicification of slot 'S' in 'fdata' ",
        "object. Indicators must be positive integers ",
        "or NA.",
        sep = ""
      ))
    }
  }
}

#' Checks consistency of slot `S` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `S` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called from the setter of slot `S` of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [setS()] for the constructor and the conditions for slots. 
".valid.S.Fdata" <- function(obj) {
  if (!all(is.na(obj@S))) {
    if (!is.integer(obj@S)) {
      stop(paste("Wrong type of slot 'S' in 'fdata' object. ",
        "Indicators must be of type 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (!all(is.na(obj@y))) {
        if (nrow(obj@S) != nrow(obj@y)) {
          stop(paste("Dimension of slot 'S' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is TRUE.",
            sep = ""
          ))
        }
      }
      if (ncol(obj@S) > 1) {
        stop(paste("Wrong dimension of slot 'S' if ",
          "slot 'bycolumn' is TRUE. Indicators ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    } else {
      if (!all(is.na(obj@y))) {
        if (ncol(obj@S) != ncol(obj@y)) {
          stop(paste("Dimension of slot 'S' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is FALSE.",
            sep = ""
          ))
        }
      }
      if (nrow(obj@S) > 1) {
        stop(paste("Wrong dimension of slot 'S' if ",
          "slot 'bycolumn' is FALSE. Indicators ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    }
    if (any(obj@S < 1)) {
      stop(paste("Wrong speicification of slot 'S' in 'fdata' ",
        "object. Indicators must be positive integers ",
        "or NA.",
        sep = ""
      ))
    }
  }
}

### Valid exp: Exposures in @exp must be of type 'numeric' or 'integer'.
### Furthermore dimensions must be conform with dimensions of data in @y.
### Exposures can only be one-dimensional and must be positive. If not
### an error is thrown.
#' Checks consistency of slot `exp` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `exp` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called during initialization of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [fdata()] for the constructor and the conditions for slots. 
".init.valid.exp.Fdata" <- function(obj) {
  if (!all(is.na(obj@exp))) {
    if (!is.numeric(obj@exp) && !is.integer(obj@exp)) {
      stop(paste("Wrong type of slot 'exp' in 'fdata' object. ",
        "Exposures must be of type 'numeric' or 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (!all(is.na(obj@y))) {
        if (nrow(obj@exp) != nrow(obj@y)) {
          warning(paste("Dimension of slot 'exp' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is TRUE",
            sep = ""
          ))
        }
      }
      if (ncol(obj@exp) > 1) {
        warning(paste("Wrong dimension of slot 'exp' if ",
          "slot 'bycolumn' is TRUE. Exposures ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    } else {
      if (!all(is.na(obj@y))) {
        if (ncol(obj@exp) != ncol(obj@y)) {
          warning(paste("Dimension of slot 'exp' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is FALSE",
            sep = ""
          ))
        }
      }
      if (nrow(obj@exp) > 1) {
        warning(paste("Wrong dimension of slot 'exp' if ",
          "slot 'bycolumn' is FALSE. Exposures ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    }
    if (any(obj@exp <= 0)) {
      stop(paste("Wrong specification of slot 'exp'. Exposures ",
        "must be positive or NA.",
        sep = ""
      ))
    }
  }
}

#' Checks consistency of slot `exp` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `exp` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called from the setter of slot `exp` of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [setExp()] for the constructor and the conditions for slots. 
".valid.exp.Fdata" <- function(obj) {
  if (!all(is.na(obj@exp))) {
    if (!is.numeric(obj@exp)) {
      stop(paste("Wrong type of slot 'exp' in 'fdata' object. ",
        "Exposures must be of type 'numeric' or 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (!all(is.na(obj@y))) {
        if (nrow(obj@exp) != nrow(obj@y)) {
          stop(paste("Dimension of slot 'exp' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is TRUE.",
            sep = ""
          ))
        }
      }
      if (ncol(obj@exp) > 1) {
        stop(paste("Wrong dimension of slot 'exp' if ",
          "slot 'bycolumn' is TRUE. Exposures ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    } else {
      if (!all(is.na(obj@y))) {
        if (ncol(obj@exp) != ncol(obj@y)) {
          stop(paste("Dimension of slot 'exp' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is FALSE.",
            sep = ""
          ))
        }
      }
      if (nrow(obj@exp) > 1) {
        stop(paste("Wrong dimension of slot 'exp' if ",
          "slot 'bycolumn' is FALSE. Exposures ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    }
    if (any(obj@exp <= 0)) {
      stop(paste("Wrong specification of slot 'exp'. Exposures ",
        "must be positive or NA.",
        sep = ""
      ))
    }
  }
}

### Valid T: Repetitions in @T must be of type 'integer'. Further,
### dimensions of @T must be consistent with @y in regard to the
### ordering in @bycolumn.
### If any element in @T is smaller than one, an error is thrown.
#' Checks consistency of slot `T` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `T` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called during initialization of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [fdata()] for the constructor and the conditions for slots. 
".init.valid.T.Fdata" <- function(obj) {
  if (!all(is.na(obj@T))) {
    if (!is.integer(obj@T)) {
      stop(paste("Wrong type of slot 'T' in 'fdata' object. ",
        "Repetitions must be of type 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (!all(is.na(obj@y))) {
        if (nrow(obj@T) != nrow(obj@y)) {
          warning(paste("Dimension of slot 'T' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is TRUE.",
            sep = ""
          ))
        }
      }
      if (ncol(obj@T) > 1) {
        warning(paste("Wrong dimension of slot 'T' if ",
          "slot 'bycolumn' is TRUE. Repetitions ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    } else {
      if (!all(is.na(obj@y))) {
        if (ncol(obj@T) != ncol(obj@y)) {
          warning(paste("Dimension of slot 'T' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is FALSE.",
            sep = ""
          ))
        }
      }
      if (nrow(obj@T) > 1) {
        warning(paste("Wrong dimension of slot 'T' if ",
          "slot 'bycolumn' is FALSE. Repetitions ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    }
    if (any(obj@T < 1)) {
      stop(paste("Wrong specification of slot 'T'. Repetitions ",
        "must be positive integers or NA.",
        sep = ""
      ))
    }
  }
}

#' Checks consistency of slot `T` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `T` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called from the setter of slot `T` of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [setT()] for the constructor and the conditions for slots. 
".valid.T.Fdata" <- function(obj) {
  if (!all(is.na(obj@T))) {
    if (!is.integer(obj@T)) {
      stop(paste("Wrong type of slot 'T' in 'fdata' object. ",
        "Repetitions must be of type 'integer'.",
        sep = ""
      ))
    }
    if (obj@bycolumn) {
      if (!all(is.na(obj@y))) {
        if (nrow(obj@T) != nrow(obj@y)) {
          stop(paste("Dimension of slot 'T' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is TRUE.",
            sep = ""
          ))
        }
      }
      if (ncol(obj@T) > 1) {
        stop(paste("Wrong dimension of slot 'T' if ",
          "slot 'bycolumn' is TRUE. Repetitions ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    } else {
      if (!all(is.na(obj@y))) {
        if (ncol(obj@T) != ncol(obj@y)) {
          stop(paste("Dimension of slot 'T' does not ",
            "match dimension of slot 'y' in 'fdata' ",
            "object if slot 'bycolumn' is FALSE.",
            sep = ""
          ))
        }
      }
      if (nrow(obj@T) > 1) {
        stop(paste("Wrong dimension of slot 'T' if ",
          "slot 'bycolumn' is FALSE. Repetitions ",
          "can only be one-dimensional.",
          sep = ""
        ))
      }
    }
    if (any(obj@T < 1)) {
      stop(paste("Wrong specification of slot 'T'. Repetitions ",
        "must be positive integers or NA.",
        sep = ""
      ))
    }
  }
}

### Valid type: The description of data type in @type must be either
### 'discrete' or 'continuous' with no exclusion. Any other choice
### throws an error.
#' Checks consistency of slot `type` of an `fdata` object
#' 
#' @description
#' Checks the consistency of slot `type` of an `fdata` object when the constructor 
#' or a setter is called. The different slots are strongly related and it has to
#'  be checked if the setting of one slot does not interfere with the definition 
#'  of another one. Is called from the setter of slot `type` of an `fdata` object.
#' 
#' @param obj An `fdata` object to be checked.
#' @returns None. Throws an error, if a certain condition is not true. 
#' @noRd
#' 
#' @seealso 
#' * [fdata] for all slots and setters of the `fdata` class.
#' * [setType()] for the constructor and the conditions for slots. 
".valid.type.Fdata" <- function(obj) {
  if (!(obj@type %in% c("discrete", "continuous"))) {
    stop(paste("Wrong choice for slot 'type'. Data can be only ",
      "'discrete' or 'continuous'",
      sep = ""
    ))
  }
}

### Valid r: The dimension of the data has to be one for 'discrete' data
### and can be greater one for 'continuous' data.
# ".valid.r.Fdata" <- function(obj)
# {
#    if (obj@type == "discrete" && obj@r > 1) {
#        stop(paste("Wrong specification of slot 'type' or slot 'r' in ",
#                   "'fdata' object. 'discrete' data can only be one-dimensional",
#                   sep = ""))
#    }
# }
