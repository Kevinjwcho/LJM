#' Create landmarked longitudinal and survival datasets for LJM/LJM
#'
#' Splits an individual-level longitudinal dataset into:
#' \itemize{
#'   \item \code{LMM_dat}: the full longitudinal record for each retained subject (all
#'   observations, including post-landmark rows). Subjects are retained if they are at
#'   risk at the landmark (\code{EvTime > s}) and have at least one longitudinal
#'   observation at or before \eqn{s} (so a landmark baseline can be defined); and
#'   \item \code{Surv_dat}: subject-level survival dataset starting at the landmark time,
#'   with \code{start = 0}, \code{stop = EvTime - s}, and event indicator carried over.
#' }
#'
#' The function also attaches \code{start/stop/event} back to \code{LMM_dat} (one row per
#' longitudinal observation) to match older LJM-style workflows.
#'
#' @param data A data frame containing longitudinal rows with at least the columns in
#'   \code{var_list} (\code{id}, \code{time}, \code{EvTime}, \code{event}).
#' @param s Numeric scalar landmark time \eqn{s}. Subjects are kept only if
#'   \code{EvTime > s}.
#' @param negative Logical. If \code{TRUE}, time is shifted to \code{time - s} and the
#'   longitudinal window is \code{[-s, 0]}. If \code{FALSE} (default), time is kept on the
#'   original scale and the window is \code{[0, s]}.
#' @param var_list A named list specifying column names in \code{data}. Defaults to
#'   \code{list(id = "id", time = "time", EvTime = "EvTime", event = "event")}.
#' @param h Numeric scalar bandwidth (half-window width) used together with
#'   \code{min_window_obs} to define the complete-case filter. The filter is
#'   applied only when both \code{h} and \code{y_vars} are supplied. Default
#'   \code{NULL}.
#' @param min_subjects Integer; the smallest number of subjects the complete-case
#'   filter may leave before the bandwidth is declared unusable. Falling below it
#'   raises a condition of class \code{"jel_bandwidth_too_small"}, which callers
#'   (e.g. bandwidth selection) can catch to skip that \code{h} rather than
#'   treating it as a generic error. Default 20.
#' @param min_window_obs Integer threshold for the complete-case filter, applied
#'   \emph{by default}: a subject is retained only if it has at least
#'   \code{min_window_obs} non-missing observations of \emph{every} marker in
#'   \code{y_vars} inside the kernel window \eqn{|t - s| \le h} (for
#'   \code{negative = TRUE} the window is centred at 0). The filter only takes
#'   effect when both \code{h} and \code{y_vars} are supplied; otherwise all
#'   at-risk subjects with an observation up to \eqn{s} are kept. Set to
#'   \code{NULL} to disable the filter even when \code{h}/\code{y_vars} are given.
#'   Default \code{2}.
#' @param y_vars Character vector of longitudinal outcome column names checked by
#'   the complete-case filter. The filter is skipped when this is \code{NULL}.
#'   Default \code{NULL}.
#'
#' @details
#' Steps performed:
#' \enumerate{
#'   \item Coerces \code{id} column to numeric.
#'   \item Keeps only subjects at risk at the landmark: \code{EvTime > s}.
#'   \item Keeps only subjects with at least one longitudinal observation at or before \eqn{s}.
#'   \item \code{LMM_dat} retains ALL longitudinal rows for those subjects (no \eqn{t \le s} cutoff).
#'   \item Builds \code{Surv_dat} with \code{start = 0}, \code{stop = EvTime - s}, and
#'   \code{event} from the subject-level record. Baseline covariates at the landmark are
#'   taken as the last observed row in the longitudinal window (\code{dplyr::last}).
#'   \item Left-joins \code{start/stop/event} back onto \code{LMM_dat}.
#' }
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{LMM_dat}: full longitudinal data for retained subjects, with \code{start/stop/event} attached.
#'   \item \code{Surv_dat}: subject-level survival data with \code{start/stop/event} and baseline covariates.
#'   \item \code{s}: the landmark time used.
#'   \item \code{var_list}: the column mapping used.
#' }
#'
#' @examples
#' \dontrun{
#' dat_split <- LJM_dat(
#'   data = dat,
#'   s = 2,
#'   negative = FALSE,
#'   var_list = list(id="id", time="time", EvTime="Time", event="event")
#' )
#' head(dat_split$LMM_dat)
#' head(dat_split$Surv_dat)
#' }
#'
#' @seealso \code{\link{LJM}}, which consumes the returned dataset.
#'
#' @export
#'
#' @importFrom dplyr group_by filter n ungroup summarise across everything last transmute left_join select any_of
#' @importFrom rlang .data
LJM_dat <- function(data, s, negative = FALSE,
                    var_list = list(id = "id", time = "time", EvTime = "EvTime", event = "event"),
                    h = NULL, min_window_obs = 2, y_vars = NULL,
                    min_subjects = 20L) {
  # min_window_obs : if set (e.g. 2), additionally keep only subjects that have at
  #   least this many observations inside the kernel window |time - s| <= h (non-NA
  #   for every outcome in y_vars).  The LJM kernel local-linear fit needs >= 2
  #   in-window points to be defined; filtering here (rather than imputing later)
  #   keeps the whole train_dataset consistent.  Requires h and y_vars.

  
  id    <- var_list[["id"]]
  time  <- var_list[["time"]]
  EvTime <- var_list[["EvTime"]]
  event <- var_list[["event"]]
  
  # id numeric
  data[[id]] <- as.numeric(data[[id]])
  
  # ----------------------------
  # Keep only subjects who are at risk at landmark s
  #   i.e., EvTime > s
  # ----------------------------
  data <- data[data[[EvTime]] > s, , drop = FALSE]

  # ----------------------------
  # Time shift if negative
  # ----------------------------
  if (negative) {
    data[[time]] <- data[[time]] - s
  }

  # ----------------------------
  # Require >= 1 longitudinal observation at or before s so that a baseline
  # covariate row at the landmark can be defined (last row at/before s).
  # data_before_s is also used to build Surv_dat baseline covariates.
  # ----------------------------
  s_cutoff <- if (negative) 0 else s

  data_before_s <- data[data[[time]] <= s_cutoff, , drop = FALSE]
  keep_ids <- unique(data_before_s[[id]])

  # ----------------------------
  # Optional >= min_window_obs in-window filter (kernel local-linear needs
  # enough local points).  Keep subjects with >= min_window_obs observations in
  # |time - center| <= h, non-NA for every outcome in y_vars.
  # ----------------------------
  if (!is.null(min_window_obs) && !is.null(h) && !is.null(y_vars)) {
    center <- if (negative) 0 else s
    dw     <- data[data[[id]] %in% keep_ids, , drop = FALSE]
    inwin  <- abs(dw[[time]] - center) <= h
    ok_ids <- keep_ids
    for (yv in y_vars) {
      cnt  <- tapply(as.integer(inwin & !is.na(dw[[yv]])), dw[[id]], sum)
      good <- as.numeric(names(cnt))[!is.na(cnt) & cnt >= min_window_obs]
      ok_ids <- intersect(ok_ids, good)
    }
    n_before <- length(keep_ids)
    keep_ids <- ok_ids
    n_after  <- length(keep_ids)

    # A bandwidth that is small relative to the visit spacing starves the local
    # fit: subjects without min_window_obs measurements inside |t - s| <= h are
    # dropped, and once too few remain the EM fails with an opaque downstream
    # error. Warn as soon as the filter bites hard, and raise a *classed*
    # condition when the retained set is unusable so that callers -- bandwidth
    # selection in particular -- can skip this h instead of treating it as a
    # generic failure.
    drop_frac <- if (n_before > 0L) 1 - n_after / n_before else 0
    if (n_after < min_subjects) {
      stop(structure(
        class = c("jel_bandwidth_too_small", "error", "condition"),
        list(message = sprintf(
          paste0("LJM_dat: bandwidth h = %g is too small at landmark s = %g: only %d of %d ",
                 "subjects have >= %d in-window observations of every marker (min_subjects = %d). ",
                 "Skip this bandwidth or widen it."),
          h, s, n_after, n_before, min_window_obs, min_subjects),
          call = NULL)))
    }
    if (drop_frac >= 0.5)
      warning(sprintf(
        paste0("LJM_dat: bandwidth h = %g at landmark s = %g drops %.0f%% of subjects ",
               "(%d of %d retained) through the min_window_obs = %d filter; ",
               "the local fit may be unstable."),
        h, s, 100 * drop_frac, n_after, n_before, min_window_obs), call. = FALSE)

    data_before_s <- data_before_s[data_before_s[[id]] %in% keep_ids, , drop = FALSE]
  }

  # ----------------------------
  # LMM_dat: longitudinal data for kept subjects (use ALL longitudinal data,
  # not just observations at or before landmark time s)
  # ----------------------------
  LMM_dat <- data[data[[id]] %in% keep_ids, , drop = FALSE]

  # ----------------------------
  # Build Surv_dat at subject level
  # - landmark time is "0"
  # - stop = EvTime - s  (time since landmark)
  # - event = event at EvTime
  # - attach baseline covariates at landmark: last row BEFORE s
  # ----------------------------
  Surv_subdat <- data_before_s[data_before_s[[id]] %in% keep_ids, , drop = FALSE] %>%
    dplyr::group_by(.data[[id]]) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::last), .groups = "drop")
  
  Surv_dat <- Surv_subdat %>%
    dplyr::transmute(
      !!id := .data[[id]],
      start = 0,
      stop  = .data[[EvTime]] - s,
      event = .data[[event]]
    ) %>%
    dplyr::left_join(
      Surv_subdat %>%
        dplyr::select(
          -dplyr::any_of(c("start", "stop", event))
        ),
      by = id
    )
  
  # ----------------------------
  # Also attach start/stop/event to LMM_dat (optional, keeps your old structure)
  # ----------------------------
  # LMM_dat <- LMM_dat %>%
  #   dplyr::left_join(Surv_dat %>% dplyr::select(.data[[id]], start, stop, event), by = id)
  # 
  LMM_dat <- LMM_dat %>%
    dplyr::select(-dplyr::any_of(c("start", "stop", "event"))) %>%  # 있으면 제거, 없으면 무시
    dplyr::left_join(
      Surv_dat %>% dplyr::select(dplyr::all_of(id), start, stop, event),
      by = id
    )
  return(list(
    LMM_dat  = LMM_dat,
    Surv_dat = Surv_dat,
    s = s,
    var_list = var_list
  ))
}