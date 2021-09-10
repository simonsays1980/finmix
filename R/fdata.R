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

### Has
### The 'hasSlot()' methods check, if the slot is not NA and returns
### TRUE if it is not NA and FALSE if it is NA.
### If argument 'verbose' is set to TRUE, an error is thrown, if
### the 'fdata' object has not the questioned slot filled.
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
### @S, @exp and @T either as column-ordered or ro-ordered matrix.
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
setMethod(
  "getY", "fdata",
  function(object) {
    return(object@y)
  }
)

setMethod(
  "getN", "fdata",
  function(object) {
    return(object@N)
  }
)

setMethod(
  "getR", "fdata",
  function(object) {
    return(object@r)
  }
)

setMethod(
  "getS", "fdata",
  function(object) {
    return(object@S)
  }
)

setMethod(
  "getBycolumn", "fdata",
  function(object) {
    return(object@bycolumn)
  }
)

setMethod(
  "getName", "fdata",
  function(object) {
    return(object@name)
  }
)

setMethod(
  "getType", "fdata",
  function(object) {
    return(object@type)
  }
)

setMethod(
  "getSim", "fdata",
  function(object) {
    return(object@sim)
  }
)

setMethod(
  "getExp", "fdata",
  function(object) {
    return(object@exp)
  }
)

setMethod(
  "getT", "fdata",
  function(object) {
    return(object@T)
  }
)

## Setters ##
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

setReplaceMethod(
  "setN", "fdata",
  function(object, value) {
    object@N <- as.integer(value)
    init.valid.Fdata(object)
    return(object)
  }
)

setReplaceMethod(
  "setR", "fdata",
  function(object, value) {
    object@r <- as.integer(value)
    .init.valid.Fdata(object)
    return(object)
  }
)

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

setReplaceMethod(
  "setName", "fdata",
  function(object, value) {
    object@name <- as.character(value)
    return(object)
  }
)

setReplaceMethod(
  "setType", "fdata",
  function(object, value) {
    object@type <- as.character(value)
    .valid.type.Fdata(object)
    return(object)
  }
)

setReplaceMethod(
  "setSim", "fdata",
  function(object, value) {
    .check.setSim.Fdata(value)
    object@sim <- value
    return(object)
  }
)

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
".plot.discrete.Fdata" <- function(obj) {
  if (has.Y(obj, verbose = TRUE)) {
    datam <- getColY(obj)
  }
  if (has.Exp(obj)) {
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
".plot.continuous.Fdata" <- function(obj, dev) {
  datam <- getColY(obj)
  if (obj@r == 1) {
    .symmetric.Hist(datam, colnames(datam))
  } else if (x@r == 2) { ## 2-dimensional
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
      expand = 0.5, lphi = 190, ltheta = 90,
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
