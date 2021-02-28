# nolint start
source("renv/activate.R")

a <- try(suppressWarnings(source(file.path(
  Sys.getenv(if (.Platform$OS.type == "windows") {
    "USERPROFILE"
  } else {
    "HOME"
  }),
  ".vscode-R",
  "init.R"
))),
silent = TRUE
)

if (interactive()) {
  suppressMessages(suppressWarnings(require(testthat)))
  suppressMessages(require(devtools))
  suppressMessages(require(usethis))
  suppressMessages(require(conflicted))
  # suppressMessages(prettycode::prettycode())

  options(
    warnPartialMatchArgs = FALSE,
    warnPartialMatchDollar = FALSE,
    warnPartialMatchAttr = FALSE,
    usethis.protocol = "https"
    # error = recover
  )


  suppressMessages(if (requireNamespace("prompt")) {
    prompt::set_prompt(function(...) {
      paste0(
        "[",
        prompt::git_branch(),
        prompt::git_dirty(),
        prompt::git_arrows(),
        "] ",
        prompt::prompt_runtime()
      )
    })
  })
}
# nolint end
