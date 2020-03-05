#' launches the shinyApp
#'
#' @export shinyAluExplorer
#'
#' @return shiny application object
#'
#' @example \dontrun {shinyAluExplorer()}
#'
#' @import shiny
#'

# wrapper for shiny::shinyApp()
shinyAluExplorer <- function() {
	# shinyAppDir("./"
	#   #system.file("shinyApp", package="meffilEWASviewer")
	# )
	shinyApp(ui = ui, server = server)
}

# NB shinyAppDir bug means globals.r is not sourced to having to use single app.r file for now.
