#' Create a new directory with a clean RMarkdown file
#'
#' This function creates a new subdirectory inside the current directory, which will
#' contain a ready-to-use RMarkdown file to be rendered in the given format.
#'
#' @param fname name of the Rmd file to create
#' @export
copy_rmd_template <- 
  function(fname = "new-doc") {
    tmp_fname <- paste0(getwd(), "/", fname, ".Rmd")
    if (file.exists(tmp_fname)) {
      stop("Cannot run copy_rmd_template()")
    }
    file.copy(
      system.file(file.path("template", "template.Rmd"), package = "ArabiDB"),
      file.path(tmp_fname)
    )
  }