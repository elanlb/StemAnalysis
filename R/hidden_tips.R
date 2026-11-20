# Estimate tips from measured, under-bark ring diameters. Columns run from
# Dnobark0 (the present profile) to the all-zero, age-zero profile.
.stem_tip_heights <- function(heights, diameters,
                              method = c("local_slope", "ratio"),
                              neighbors = 3L) {
  method <- match.arg(method)
  if (length(neighbors) != 1L || !is.numeric(neighbors) ||
      !is.finite(neighbors) || neighbors < 1 || neighbors != floor(neighbors)) {
    stop("tip_neighbors must be a positive integer")
  }
  if (!is.numeric(heights) || length(heights) < 2L ||
      any(!is.finite(heights)) || any(diff(heights) <= 0)) {
    stop("Disc heights must be finite and strictly increasing")
  }
  diameters <- as.matrix(diameters)
  if (!is.numeric(diameters) || nrow(diameters) != length(heights) ||
      ncol(diameters) < 3L || any(!is.finite(diameters)) || any(diameters < 0)) {
    stop("Ring diameters must be a finite, non-negative numeric matrix")
  }
  if (any(diameters[nrow(diameters), ] != 0) ||
      any(diameters[, ncol(diameters)] != 0)) {
    stop("Include an all-zero top row and an all-zero age-zero column")
  }
  if (any(diameters[, -1L, drop = FALSE] >
          diameters[, -ncol(diameters), drop = FALSE])) {
    stop("Ring diameters must be nested from Dnobark0 to the oldest ring")
  }
  for (j in seq_len(ncol(diameters) - 1L)) {
    if (any(diff(as.integer(diameters[, j] > 0)) > 0)) {
      stop("A ring cannot reappear above a disc where it is absent")
    }
  }

  tips <- numeric(ncol(diameters) - 2L)
  enclosing_tip <- heights[length(heights)]
  for (k in seq_along(tips)) {
    column <- k + 1L
    upper <- which(diameters[, column] == 0)[1L]
    lower <- upper - 1L
    if (lower == 0L) {
      tips[k] <- heights[1L]
      enclosing_tip <- tips[k]
      next
    }
    interval <- heights[upper] - heights[lower]
    diameter <- diameters[lower, column]
    outer_diameter <- diameters[lower, column - 1L]

    if (method == "ratio") {
      # Retain the original height estimator for explicit comparisons.
      tips[k] <- heights[lower] + diameter / outer_diameter * interval
      next
    }

    outer <- seq_len(column - 1L)
    drop <- diameters[lower, outer] - diameters[upper, outer]
    measured <- outer[diameters[upper, outer] > 0 & drop > 0]
    if (length(measured)) {
      # Largest column indices are the closest enclosing measured rings.
      measured <- utils::tail(measured, min(neighbors, length(measured)))
      taper <- stats::median((diameters[lower, measured] -
                              diameters[upper, measured]) / interval)
      tip <- heights[lower] + diameter / taper
    } else if (upper == nrow(diameters)) {
      # No ring crosses the final section: use the observed present-day cone.
      taper <- diameters[lower, 1L] / interval
      tip <- heights[lower] + diameter / taper
    } else {
      # Flat or outward slopes cannot define an inward tip. Use a bounded
      # diameter ratio, anchored to the already estimated enclosing tip.
      top <- min(heights[upper], enclosing_tip)
      tip <- heights[lower] + diameter / outer_diameter *
        (top - heights[lower])
    }

    # Respect absence at the upper disc and keep older tips inside newer ones.
    tips[k] <- max(heights[lower],
                   min(tip, heights[upper], enclosing_tip))
    enclosing_tip <- tips[k]
  }
  tips
}
