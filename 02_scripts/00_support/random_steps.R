#' Generate Random Steps

.random_steps.bursted_steps_xyt <- function(
    x, n_control = 10,
    sl_distr = fit_distr(x$sl_, "gamma"), # this argument could be remove
    ta_distr = fit_distr(x$ta_, "vonmises"), # this argument could be remove
    rand_sl = random_numbers(sl_distr, n = 1e5),
    rand_ta = random_numbers(ta_distr, n = 1e5),
    include_observed = TRUE, ...) {
  
  bursts <- split(x, x$burst_)
  
  if (any(len.ok <- sapply(bursts, nrow) < 3)) {
    warning("Some bursts contain < 3 steps and will be removed")
    bursts <- bursts[!len.ok]
  }
  
  start_ids <- c(1, head(cumsum(rle(unlist(sapply(bursts, "[[", "burst_")))$lengths), -1) + 1)
  
  pb <- progress::progress_bar$new(total = length(bursts))
  out <- lapply(seq_along(bursts), function(i) {
    q <- bursts[[i]]
    pb$tick()
    class(q) <- class(q)[-1]
    if (nrow(q) > 1) {
      random_steps(q, n_control = n_control, sl_distr = NULL,
                   ta_distr = NULL, rand_sl = rand_sl,
                   rand_ta = rand_ta, include_observed = include_observed,
                   start_id = start_ids[i], ...)
    }
  })
  
  out <- out[!sapply(out, is.null)]
  
  cls <- class(out[[1]])
  
  out <- data.table::rbindlist(
    out
  ) |> tibble::as_tibble()
  
  class(out) <- cls
  attributes(out)$sl_ <- sl_distr
  attributes(out)$ta_ <- ta_distr
  attributes(out)$n_control_ <- n_control
  attr(out, "crs_") <- attr(x, "crs_")
  out
}



#' @export
plot.random_steps <- function(x, ...) {
  plot(0, 0, type = "n",
       xlim = grDevices::extendrange(c(x$x1_, x$x2_), f = 0.1),
       ylim = grDevices::extendrange(c(x$y1_, x$y2_), f = 0.1),
       xlab = "x", ylab = "y")
  for (i in 1:nrow(x)) {
    graphics::lines(c(x$x1_[i], x$x2_[i]), c(x$y1_[i], x$y2_[i]), lty = 2,
                    col = "grey79")
  }
  
  x1 <- x[x$case_, ]
  for (i in 1:nrow(x1)) {
    lines(c(x1$x1_[i], x1$x2_[i]), c(x1$y1_[i], x1$y2_[i]))
    graphics::points(x1$x1_[i], x1$y1_[i], pch = 20, cex = 2, col = "red")
    graphics::points(x1$x2_[i], x1$y2_[i], pch = 20, cex = 2, col = "red")
  }
}

# Flag incomplete steps ----

#' Remove strata with missing values for observed steps
#'
#' @param x An object of class `random_steps`.
#' @param col A character with the column name that will be scanned for missing values.
#' @template dots_none
#' @name remove_incomplete_strata
#'
#' @return An object of class `random_steps`, where observed steps (`case_ == TRUE`) with missing values (`NA`) in the column `col` are removed (including all random steps).

#' @examples
#'
#' mini_deer <- deer[1:4, ]
#'
#' # The first step is removed, because we have `NA` turn angles.
#' mini_deer |> steps() |> random_steps() |> remove_incomplete_strata() |>
#'   select(case_, ta_, step_id_)

#' @export
remove_incomplete_strata <- function(x, ...) {
  UseMethod("remove_incomplete_strata", x)
}

#' @export
#' @rdname remove_incomplete_strata
remove_incomplete_strata.random_steps <- function(x, col = "ta_", ...) {
  
  checkmate::assert_character(col, len = 1)
  
  if (!col %in% names(x)) {
    stop("`col` not found in `x` (make sure the column name is spelled correct).")
    
  }
  
  x.case <- dplyr::filter(x, case_)
  incomplete.steps <- which(is.na(x.case[[col]]))
  dplyr::filter(
    x, !step_id_ %in% x.case$step_id_[incomplete.steps])
}




rsteps_transfer_attr <- function(from, to) {
  from <- attributes(from)
  attributes(to)$class <- from$class
  attributes(to)$sl_ <- from$sl_
  attributes(to)$ta_ <- from$ta_
  attributes(to)$crs_ <- from$crs_
  to
}



# see here: https://github.com/hadley/dplyr/issues/719
#' @export
arrange.random_steps <- function(.data, ..., .dots) {
  xx <- NextMethod()
  rsteps_transfer_attr(.data, xx)
}


#' @export
filter.random_steps <- function(.data, ..., .dots) {
  xx <- NextMethod()
  rsteps_transfer_attr(.data, xx)
}

#' @export
group_by.random_steps <- function(.data, ..., .dots) {
  xx <- NextMethod()
  rsteps_transfer_attr(.data, xx)
}

#' @export
nest.random_steps <- function(.data, ..., .dots) {
  NextMethod()
}

#' @export
select.random_steps <- function(.data, ..., .dots) {
  xx <- NextMethod()
  rsteps_transfer_attr(.data, xx)
}

#' @export
summarise.random_steps <- function(.data, ..., .dots) {
  NextMethod()
}


#' @export
summarize.random_steps <- function(.data, ..., .dots) {
  NextMethod()
}