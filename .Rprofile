source("renv/activate.R")

if (interactive()) {
  suppressMessages(require(devtools))
  suppressMessages(require(usethis))
  # suppressMessages(require(conflicted))
  # suppressMessages(prettycode::prettycode())

  options(
    warnPartialMatchArgs = FALSE,
    warnPartialMatchDollar = FALSE,
    warnPartialMatchAttr = FALSE,
    usethis.protocol = "https"
  )

  prompt::set_prompt(function(...) {
    paste0(
      "[", prompt::git_branch(), prompt::git_dirty(), prompt::git_arrows(), "] ",
      prompt::prompt_runtime()
    )
  })
}


