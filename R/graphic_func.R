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

### Private functions.
### These functions are not exported.

### Checking
### This function checks, if an option 'title' for the
### graphical device used by R is available. If the answer
### is TRUE, the title can be set by a 'plot()' function.
".check.grDevice" <- function() {
  ## title argument ##
  any(names(formals(getOption("device")))
  == "title")
}

### Plotting
### This functions checks the dimension of a dataset 'y'
### an distributes histograms for each variable in the
### dataset symmetrically around the graphical grid.
".symmetric.Hist" <- function(y, lab.names) {
  r <- NCOL(y)
  if (r == 1) {
    .comb.Hist(y, lab.names)
  } else if (r == 2) {
    par(
      mfrow = c(1, 2), mar = c(2, 2, 2, 2),
      oma = c(4, 5, 1, 5)
    )
    for (i in 1:2) {
      .comb.Hist(y[, i], lab.names[i])
    }
  } else if (r == 3) {
    par(
      mfrow = c(1, 3), mar = c(2, 2, 2, 2),
      oma = c(4, 5, 1, 5)
    )
    for (i in 1:r) {
      .comb.Hist(y[, i], lab.names[i])
    }
  } else if (r > 3 && r < 17 && sqrt(r) %% 1 == 0) {
    par(
      mfrow = c(sqrt(r), sqrt(r)),
      mar = c(2, 2, 2, 2),
      oma = c(4, 5, 1, 5)
    )
    for (i in 1:r) {
      .comb.Hist(y[, i], lab.names[i])
    }
  } else {
    if (r == 5) {
      par(
        mfrow = c(2, 3),
        mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
      for (i in 1:4) {
        .comb.Hist(y[, i], lab.names[i])
      }
      plot.new()
      .comb.Hist(y[, r], lab.names[r])
    } else if (r == 6) {
      par(
        mfrow = c(2, 3),
        mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
      for (i in 1:r) {
        .comb.Hist(y[, i], lab.names[i])
      }
    } else {
      if (r %% 2 == 0) {
        ## check how many rows can be completely
        ## filled
        n <- r %/% 4
        par(
          mfrow = c(n, 4),
          mar = c(2, 2, 2, 2),
          oma = c(4, 5, 1, 5)
        )
        for (i in 1:(n * 4)) {
          .comb.Hist(y[, i], lab.names[i])
        }
        ## if some rows cannot be completely
        ## filled, fill them symmetrically
        ## there can only be two left:
        .comb.Hist(y[, r - 1], lab.names[r - 1])
        plot.new()
        .comb.Hist(y[, r], lab.names[r])
      } else {
        n <- r %/% 5
        par(
          mfrow = c(n, 5),
          mar = c(2, 2, 2, 2),
          oma = c(4, 5, 1, 5)
        )
        for (i in 1:(n * 5)) {
          .comb.Hist(y[, i], lab.names[i])
        }
        ## if some rows cannot be completely,
        ## filled, fill them symmetrically
        ## either there are two left or four
        ## left
        if (r %% 5 == 2) {
          plot.new()
          .comb.Hist(y[, r - 1], lab.names[r - 1])
          plot.new()
          .comb.Hist(y[, r], lab.names[r])
          plot.new()
        } else {
          .comb.Hist(y[, r - 3], lab.names[r - 3])
          .comb.Hist(y[, r - 2], lab.names[r - 2])
          plot.new()
          .comb.Hist(y[, r - 1], lab.names[r - 1])
          .comb.Hist(y[, r], lab.names[r])
        }
      }
    }
  }
}

### This functions checks the dimension of a dataset 'y'
### an distributes Kernel densities for each variable in the
### dataset symmetrically around the graphical grid.
".symmetric.Dens" <- function(y, lab.names) {
  r <- NCOL(y)
  if (r == 1) {
    .comb.Dens(y, lab.names)
  } else if (r == 2) {
    par(
      mfrow = c(1, 2), mar = c(2, 2, 2, 2),
      oma = c(4, 5, 1, 5)
    )
    for (i in 1:2) {
      .comb.Dens(y[, i], lab.names[i])
    }
  } else if (r == 3) {
    par(
      mfrow = c(1, 3), mar = c(2, 2, 2, 2),
      oma = c(4, 5, 1, 5)
    )
    for (i in 1:r) {
      .comb.Dens(y[, i], lab.names[i])
    }
  } else if (r > 3 && r < 17 && sqrt(r) %% 1 == 0) {
    par(
      mfrow = c(sqrt(r), sqrt(r)),
      mar = c(2, 2, 2, 2),
      oma = c(4, 5, 1, 5)
    )
    for (i in 1:r) {
      .comb.Dens(y[, i], lab.names[i])
    }
  } else {
    if (r == 5) {
      par(
        mfrow = c(2, 3),
        mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
      for (i in 1:4) {
        .comb.Dens(y[, i], lab.names[i])
      }
      plot.new()
      .comb.Dens(y[, r], lab.names[r])
    } else if (r == 6) {
      par(
        mfrow = c(2, 3),
        mar = c(2, 2, 2, 2),
        oma = c(4, 5, 1, 5)
      )
      for (i in 1:r) {
        .comb.Dens(y[, i], lab.names[i])
      }
    } else {
      if (r %% 2 == 0) {
        ## check how many rows can be completely
        ## filled
        n <- r %/% 4
        par(
          mfrow = c(n, 4),
          mar = c(2, 2, 2, 2),
          oma = c(4, 5, 1, 5)
        )
        for (i in 1:(n * 4)) {
          .comb.Dens(y[, i], lab.names[i])
        }
        ## if some rows cannot be completely
        ## filled, fill them symmetrically
        ## there can only be two left:
        .comb.Dens(y[, r - 1], lab.names[r - 1])
        plot.new()
        .comb.Dens(y[, r], lab.names[r])
      } else {
        n <- r %/% 5
        par(
          mfrow = c(n, 5),
          mar = c(2, 2, 2, 2),
          oma = c(4, 5, 1, 5)
        )
        for (i in 1:(n * 5)) {
          .comb.Dens(y[, i], lab.names[i])
        }
        ## if some rows cannot be completely,
        ## filled, fill them symmetrically
        ## either there are two left or four
        ## left
        if (r %% 5 == 2) {
          plot.new()
          .comb.Dens(y[, r - 1], lab.names[r - 1])
          plot.new()
          .comb.Dens(y[, r], lab.names[r])
          plot.new()
        } else {
          .comb.Dens(y[, r - 3], lab.names[r - 3])
          .comb.Dens(y[, r - 2], lab.names[r - 2])
          plot.new()
          .comb.Dens(y[, r - 1], lab.names[r - 1])
          .comb.Dens(y[, r], lab.names[r])
        }
      }
    }
  }
}

### This function plots a histogram with 'finmix' specific
### settings. In addition it uses 'rug()' to plot the data
### points.
".comb.Hist" <- function(y, lab.name) {
  hist(y,
    col = "gray65",
    border = "white", cex = 0.7,
    cex.axis = 0.7, freq = TRUE,
    xlab = "", main = "", cex.lab = 0.7
  )
  rug(y, col = "gray47")
  mtext(
    side = 1, do.call(bquote, lab.name),
    cex = 0.7, line = 3
  )
}

### This function plots a Kernel density with 'finmix' specific
### settings. In addition it uses 'rug()' to plot the data
### points.
".comb.Dens" <- function(y, lab.name) {
  dens <- bkde(y)
  plot(dens$x, dens$y,
    col = "gray47",
    cex.axis = .7, cex = .7, type = "l",
    xlab = "", main = "", ylab = "Density",
    cex.lab = .7
  )
  rug(y, col = "gray47")
  mtext(
    side = 1, do.call(bquote, lab.name),
    cex = 0.7, line = 3
  )
}
