# Part of nlmixr2sir.
# Run identity: fingerprinting, directory ownership, and the on-disk manifest.

# Bump when the shape of the saved state changes incompatibly. Recovery refuses
# a state written by a different version rather than guessing at its layout.
.sirStateVersion <- 2L

# Bump when the numerical algorithm changes in a way that makes an older saved
# run no longer an answer to the same question. Distinct from the state version,
# which is only about the on-disk layout, and from the package version, which
# moves for documentation changes too.
.sirAlgorithmVersion <- 1L

.sirStateSchema <- function() {
  list(prefix = "sir", version = .sirStateVersion)
}

.sirMarkerFile <- "sir_manifest.dcf"

# Controls that change the answer. Parallelism and resume behaviour do not, so
# a run must not be invalidated by having been resumed with different worker
# settings -- the seeding is per iteration and independent of them.
.sirStatisticalControls <- c(
  "thetaInflation", "omegaInflation", "sigmaInflation",
  "capCorrelation", "capResampling", "recenter", "boxcox",
  "omegaFallback", "sigmaFallbackRse", "omegaDf",
  "rseTheta", "rseOmega", "rseSigma",
  "covmatInput", "rawresInput", "offsetRawres", "inFilter",
  "objfTolerance"
)

# A short, stable digest of an arbitrary R object.
#
# Base R has no object hash, so the object is serialized to a temporary file and
# md5-summed. Version 3 serialization of plain data is deterministic, which is
# all that is needed here: the digest only ever has to agree with itself.
# Anything that cannot be serialized digests to NA, and an NA field is treated
# as "cannot verify" rather than as a mismatch.
.sirDigest <- function(x) {
  tryCatch(
    {
      tmp <- tempfile("sirfp")
      on.exit(unlink(tmp), add = TRUE)
      con <- file(tmp, open = "wb")
      serialize(x, con, version = 3L, ascii = FALSE)
      close(con)
      unname(tools::md5sum(tmp))
    },
    error = function(e) NA_character_
  )
}

# Everything that, if changed, means a saved run is no longer the run being
# asked for. Stored alongside the state so recovery can refuse a mismatch and
# name the field that moved, rather than returning a stale result under a new
# label.
# A schedule frame from the two vectors runSIR() is called with. One
# constructor so the result attribute, the saved state and the fingerprint
# cannot drift into different shapes.
.sirSchedule <- function(nSamples, nResample) {
  data.frame(
    iter = seq_along(nSamples),
    nSamples = as.integer(nSamples),
    nResample = as.integer(nResample)
  )
}

# `schedule` is a data frame of iter/nSamples/nResample covering every iteration
# the run will contain -- for an extended run that is the prior schedule plus
# the extension, not the extension alone. Taking the frame rather than the two
# vectors is what keeps the saved identity and the result's schedule attribute
# from drifting apart.
.sirRunFingerprint <- function(fit, ps, schedule, control, initial = NULL) {
  # fit$ui$funTxt is the model block as text, which is what distinguishes one
  # structural model from another. fit$uiFun is NULL on a fitted object, so
  # deparsing it yields the string "NULL" for every fit -- which silently made
  # every model look identical.
  model_src <- tryCatch(
    paste(as.character(fit$ui$funTxt), collapse = "\n"),
    error = function(e) NA_character_
  )
  if (length(model_src) != 1L || is.na(model_src) || !nzchar(model_src)) {
    model_src <- tryCatch(
      paste(deparse(fit$ui$fun), collapse = "\n"),
      error = function(e) NA_character_
    )
  }
  if (
    length(model_src) != 1L ||
      is.na(model_src) ||
      !nzchar(model_src) ||
      identical(model_src, "NULL")
  ) {
    model_src <- NA_character_
  }

  dat <- tryCatch(nlmixr2est::getData(fit), error = function(e) NULL)
  data_rows <- if (is.null(dat)) NA_integer_ else nrow(dat)
  data_digest <- if (is.null(dat)) {
    NA_character_
  } else {
    # Strip attributes that vary without the data varying: the class carried by
    # different table backends, row names, and so on.
    .sirDigest(lapply(as.data.frame(dat), as.vector))
  }

  # The numbers the run actually proposes from. Digesting the RESOLVED proposal
  # covers every input route at once -- fit$cov, the parsed contents of a
  # covmatInput file, and the selected rows of a rawresInput file -- where
  # digesting the control object only ever captured a path string, so replacing
  # the file at that path left the fingerprint unchanged.
  proposal_digest <- if (is.null(initial)) {
    NA_character_
  } else {
    .sirDigest(list(
      cov = round(unname(as.matrix(initial$covMat)), 12L),
      mu = round(unname(as.numeric(initial$mu %||% numeric(0))), 12L),
      source = as.character(initial$source %||% NA_character_)
    ))
  }

  # Names alone do not pin the estimation problem: the same names with
  # different kinds or bounds is a different problem.
  schema_digest <- .sirDigest(lapply(
    ps[intersect(c("sirName", "kind", "lower", "upper", "neta1", "neta2"), names(ps))],
    as.vector
  ))

  list(
    stateVersion = .sirStateVersion,
    algoVersion = .sirAlgorithmVersion,
    model = .sirDigest(model_src),
    dataRows = data_rows,
    data = data_digest,
    params = paste(ps$sirName, collapse = "|"),
    paramSchema = schema_digest,
    proposal = proposal_digest,
    estimates = .sirDigest(round(as.numeric(ps$est), 10L)),
    objf = .sirDigest(round(as.numeric(fit$objf), 8L)),
    # .sirFitEst(), not fit$est: on a fit whose tables carry an `est` column
    # the latter is a per-row vector, which would land in the fingerprint and
    # the manifest as 132 repeats of the method name.
    estMethod = .sirFitEst(fit),
    # The method the NUMBERS came from, which is not always the method the fit
    # was run with: the imp family is scored as FOCEi (see .sirEvalMethod()).
    # estMethod alone would describe such a run as "impmap" when every objective
    # in it was produced by FOCEi.
    #
    # Derived from estMethod today, so it adds no discriminating power right
    # now. It is an identity field anyway, so that a state file written under a
    # different mapping is refused rather than silently reused.
    evalMethod = .sirEvalMethod(fit),
    schedule = paste0(
      paste(as.integer(schedule$nSamples), collapse = ","),
      "/",
      paste(as.integer(schedule$nResample), collapse = ",")
    ),
    controls = .sirDigest(control[.sirStatisticalControls]),
    pkgVersions = paste(
      vapply(
        c("nlmixr2sir", "nlmixr2est", "nlmixr2utils", "rxode2"),
        function(p) {
          tryCatch(
            as.character(utils::packageVersion(p)),
            error = function(e) "?"
          )
        },
        character(1L)
      ),
      collapse = "/"
    )
  )
}

# Fields whose disagreement means the saved run is a different run. Package
# versions are recorded but not enforced: a dependency bump does not by itself
# invalidate a result, and the objective preflight catches it if it did.
.sirIdentityFields <- c(
  "stateVersion", "algoVersion", "model", "dataRows", "data", "params",
  "paramSchema", "proposal", "estimates", "objf", "estMethod", "evalMethod",
  "schedule", "controls"
)

# Returns the names of fields that disagree. `ignore` exempts fields that are
# expected to differ: addIterations deliberately changes the schedule.
# Returns the names of fields that disagree, with the subset that could not be
# verified at all carried in an "unverifiable" attribute.
#
# `ignore` exempts fields expected to differ -- addIterations deliberately
# changes the schedule.
#
# `failClosed` decides what an undigestible field means. Skipping it is the
# right answer for a fresh run: an unusual fit should not be blocked because one
# field would not serialize. It is the wrong answer for recovery, where the
# whole point is to establish that the saved run IS this run. Skipping there
# means reuse is permitted precisely when identity cannot be shown.
.sirCompareFingerprints <- function(saved, current, ignore = character(),
                                    failClosed = FALSE) {
  if (!is.list(saved)) {
    out <- "stateVersion"
    attr(out, "unverifiable") <- "stateVersion"
    return(out)
  }
  fields <- setdiff(.sirIdentityFields, ignore)
  bad <- character()
  unverifiable <- character()
  missingOk <- function(v) is.null(v) || (length(v) == 1L && is.na(v))
  for (f in fields) {
    a <- saved[[f]]
    b <- current[[f]]
    if (missingOk(a) || missingOk(b)) {
      # The state version is never optional: absent means a state this code
      # does not understand.
      if (identical(f, "stateVersion")) {
        bad <- c(bad, f)
        unverifiable <- c(unverifiable, f)
        next
      }
      unverifiable <- c(unverifiable, f)
      if (failClosed) {
        bad <- c(bad, f)
      }
      next
    }
    if (!isTRUE(all.equal(a, b))) {
      bad <- c(bad, f)
    }
  }
  attr(bad, "unverifiable") <- unverifiable
  bad
}

.sirFingerprintLabels <- c(
  stateVersion = "saved-state format version",
  algoVersion = "SIR algorithm version",
  paramSchema = "parameter kinds and bounds",
  proposal = "resolved initial proposal",
  model = "model definition",
  dataRows = "number of data rows",
  data = "dataset contents",
  params = "estimated parameter set",
  estimates = "parameter estimates",
  objf = "objective function value",
  estMethod = "estimation method",
  schedule = "sample/resample schedule",
  controls = "statistical control settings"
)

.sirAbortFingerprint <- function(bad, dir, action) {
  labels <- unname(.sirFingerprintLabels[bad])
  labels[is.na(labels)] <- bad[is.na(labels)]
  unver <- attr(bad, "unverifiable", exact = TRUE)
  unverLabels <- if (length(unver) > 0L) {
    l <- unname(.sirFingerprintLabels[unver])
    l[is.na(l)] <- unver[is.na(l)]
    l
  } else {
    character()
  }
  cli::cli_abort(c(
    "The saved SIR run in {.path {dir}} does not match this {action} request.",
    "x" = "Changed: {.val {labels}}.",
    if (length(unverLabels) > 0L) {
      c("x" = "Could not be verified at all: {.val {unverLabels}}.")
    },
    "i" = "Returning the stored result would label one run's output as another's.",
    "i" = "Use a different {.arg directory}, or {.code runSIRControl(recover = FALSE)} to start fresh."
  ))
}

# --- Directory ownership -----------------------------------------------------

# Read the manifest back, or NULL if there isn't a readable one.
.sirReadManifest <- function(dir) {
  path <- file.path(dir, .sirMarkerFile)
  if (!file.exists(path)) {
    return(NULL)
  }
  m <- tryCatch(read.dcf(path), error = function(e) NULL)
  if (is.null(m) || nrow(m) < 1L || ncol(m) < 1L) {
    return(NULL)
  }
  stats::setNames(as.list(as.character(m[1L, ])), colnames(m))
}

# Ownership is asserted by a manifest this package wrote, and it is checked by
# reading the manifest -- not by testing for the filename.
#
# The filename alone is not proof: an empty file, a file written by another
# tool, or one describing a state version this code does not understand would
# all have passed. Since the marker is what authorises recursive deletion,
# accepting any of those turns an unrelated directory into a deletable one.
.sirDirIsOwned <- function(dir) {
  man <- .sirReadManifest(dir)
  if (is.null(man)) {
    return(FALSE)
  }
  if (!identical(man[["Package"]], "nlmixr2sir")) {
    return(FALSE)
  }
  if (!identical(man[["Prefix"]], "sir")) {
    return(FALSE)
  }
  sv <- suppressWarnings(as.integer(man[["StateVersion"]]))
  if (length(sv) != 1L || is.na(sv) || sv > .sirStateVersion) {
    return(FALSE)
  }
  TRUE
}

.sirDirIsEmpty <- function(dir) {
  length(list.files(dir, all.files = TRUE, no.. = TRUE)) == 0L
}

# The run manifest doubles as the ownership marker. Human-readable on purpose:
# it records what produced the directory, and it is what later permits the
# directory to be overwritten.
.sirWriteManifest <- function(dir, fingerprint, fitName, seed = NULL) {
  fields <- list(
    Package = "nlmixr2sir",
    Version = tryCatch(
      as.character(utils::packageVersion("nlmixr2sir")),
      error = function(e) "?"
    ),
    Prefix = "sir",
    Fit = as.character(fitName %||% NA_character_),
    Created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    StateVersion = as.character(.sirStateVersion),
    Schedule = as.character(fingerprint$schedule %||% NA_character_),
    Parameters = as.character(fingerprint$params %||% NA_character_),
    EstimationMethod = as.character(fingerprint$estMethod %||% NA_character_),
    DataRows = as.character(fingerprint$dataRows %||% NA_character_),
    PackageVersions = as.character(fingerprint$pkgVersions %||% NA_character_),
    Seed = if (is.null(seed)) NA_character_ else as.character(seed)
  )
  df <- as.data.frame(fields, stringsAsFactors = FALSE, check.names = FALSE)
  # Fatal rather than a warning: the manifest is what authorises this directory
  # to be cleared later, so a run that could not write one must not proceed as
  # though it had.
  path <- file.path(dir, .sirMarkerFile)
  tryCatch(
    write.dcf(df, path),
    error = function(e) {
      cli::cli_abort(c(
        "Could not write the SIR run manifest to {.path {path}}.",
        "x" = conditionMessage(e),
        "i" = "The manifest records what produced the directory and gates any later overwrite, so the run stops here."
      ))
    }
  )
  invisible(path)
}

# May this run write into `dir` at all?
#
# Three cases are legitimate: a directory this call created, an empty one, and
# one this package already owns. Anything else is somebody's data.
#
# This has to be checked BEFORE the seed, state, or manifest is written. The
# earlier code guarded only the overwrite path, so with the default
# recover = TRUE an existing non-empty directory came back in resume mode, held
# no SIR state, and had a manifest written into it anyway -- which made it
# "owned", and a later fresh run was then free to delete everything in it.
# Writing the marker must not be what creates the ownership it checks for.
.sirAssertClaimable <- function(dir) {
  if (!dir.exists(dir) || .sirDirIsEmpty(dir) || .sirDirIsOwned(dir)) {
    return(invisible(TRUE))
  }
  marker <- .sirMarkerFile # nolint: object_usage_linter.
  cli::cli_abort(c(
    "{.path {dir}} is not empty and does not appear to be a nlmixr2sir run directory.",
    "x" = "Refusing to write run artifacts into it.",
    "i" = "A directory created by {.fn runSIR} carries a valid {.file {marker}} manifest; this one does not.",
    "i" = "Point {.arg directory} at a new or empty location."
  ))
}

# resolveRunDir() returns overwrite mode for any existing directory the user
# named explicitly. That is not enough to justify unlink(recursive = TRUE): an
# explicitly supplied directory can hold anything at all.
.sirAssertSafeToClear <- function(dir) {
  if (!dir.exists(dir) || .sirDirIsEmpty(dir) || .sirDirIsOwned(dir)) {
    return(invisible(TRUE))
  }
  # Bound locally so cli's glue can see it. # nolint: object_usage_linter.
  marker <- .sirMarkerFile
  cli::cli_abort(c(
    "{.path {dir}} is not empty and does not appear to be a nlmixr2sir run directory.",
    "x" = "Refusing to delete its contents.",
    "i" = "A directory created by {.fn runSIR} carries a {.file {marker}} manifest; this one does not.",
    "i" = "Point {.arg directory} somewhere else, or remove it yourself if you are sure."
  ))
}
