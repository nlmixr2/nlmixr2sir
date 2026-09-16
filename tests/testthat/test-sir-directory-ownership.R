# B3: a directory must be owned before it can be claimed or cleared
#
# `.sirAssertSafeToClear()` guarded the overwrite path only. With the default
# `recover = TRUE`, `resolveRunDir()` returns resume mode for any existing
# directory, `saved_state` comes back NULL when it holds no SIR state, and the
# run wrote its manifest anyway -- silently marking an unrelated directory as
# nlmixr2sir's. A later fresh run then finds a valid-looking marker and is
# authorised to delete everything in it.
#
# Ownership must be established before the marker is written, not created by
# writing it.

test_that("a non-empty unowned directory is refused, not claimed", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  sentinel <- file.path(dir, "sentinel.txt")
  writeLines("precious", sentinel)

  err <- tryCatch(
    .sirQuiet(runSIR(
      theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(recover = TRUE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "not appear to be|nlmixr2sir")

  # Nothing was written: no marker, no seed, no state.
  expect_false(.sirDirIsOwned(dir))
  expect_setequal(list.files(dir, all.files = TRUE, no.. = TRUE), "sentinel.txt")
  expect_equal(readLines(sentinel), "precious")
})

test_that("the two-call sequence cannot destroy an unrelated file", {
  skip_on_cran()
  # The path the review described: a failed recovery attempt that leaves a
  # marker behind, followed by a fresh run that clears the directory.
  dir <- withr::local_tempdir()
  sentinel <- file.path(dir, "sentinel.txt")
  writeLines("precious", sentinel)

  try(.sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(recover = TRUE, workers = 1L)
  )), silent = TRUE)

  try(.sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(recover = FALSE, workers = 1L)
  )), silent = TRUE)

  expect_true(file.exists(sentinel))
  expect_equal(readLines(sentinel), "precious")
})

test_that("an empty directory may be claimed", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0L)
  expect_no_error(.sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(recover = TRUE, workers = 1L)
  )))
  expect_true(.sirDirIsOwned(dir))
})

test_that("ownership requires a valid manifest, not just the filename", {
  skip_on_cran()
  dir <- withr::local_tempdir()

  # An empty file with the right name proves nothing.
  file.create(file.path(dir, "sir_manifest.dcf"))
  expect_false(.sirDirIsOwned(dir))

  # Nor does a well-formed DCF written by something else.
  write.dcf(
    data.frame(Package = "someone.else", Prefix = "scm", stringsAsFactors = FALSE),
    file.path(dir, "sir_manifest.dcf")
  )
  expect_false(.sirDirIsOwned(dir))

  # Nor a manifest claiming a state version this code does not understand.
  write.dcf(
    data.frame(
      Package = "nlmixr2sir", Prefix = "sir",
      StateVersion = "999", stringsAsFactors = FALSE
    ),
    file.path(dir, "sir_manifest.dcf")
  )
  expect_false(.sirDirIsOwned(dir))
})

test_that("a foreign manifest does not authorize deletion", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  sentinel <- file.path(dir, "sentinel.txt")
  writeLines("precious", sentinel)
  write.dcf(
    data.frame(Package = "someone.else", Prefix = "scm", stringsAsFactors = FALSE),
    file.path(dir, "sir_manifest.dcf")
  )

  expect_error(.sirAssertSafeToClear(dir), "not appear to be|nlmixr2sir")
  expect_true(file.exists(sentinel))
})

test_that("a manifest this package wrote does authorize deletion", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  writeLines("scratch", file.path(dir, "leftover.txt"))
  .sirWriteManifest(
    dir,
    .sirRunFingerprint(
      theoFit(), .sirParamSpace(theoFit()), .sirSchedule(16L, 8L), runSIRControl(workers = 1L)
    ),
    fitName = "theoFit"
  )
  expect_true(.sirDirIsOwned(dir))
  expect_silent(.sirAssertSafeToClear(dir))
})

test_that("a manifest that cannot be written is fatal", {
  skip_on_cran()
  # The marker gates recursive deletion, so a run that could not write one must
  # not continue as though it had.
  dir <- withr::local_tempdir()
  fp <- .sirRunFingerprint(
    theoFit(), .sirParamSpace(theoFit()), .sirSchedule(16L, 8L), runSIRControl(workers = 1L)
  )
  # suppressWarnings() covers base R's own "cannot open file" warning from
  # write.dcf, which fires on the way to the error. The error is the subject
  # here; the warning is incidental and is not a condition this package raises.
  expect_error(
    suppressWarnings(.sirWriteManifest(
      file.path(dir, "no", "such", "directory"), fp, fitName = "theoFit"
    )),
    "manifest"
  )
})
