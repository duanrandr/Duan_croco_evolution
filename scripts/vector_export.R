# Step 1: Use standard installed packages; do not bundle platform libraries.
if (!requireNamespace('svglite', quietly=TRUE)) {
  stop("Install the svglite package before building figures: install.packages('svglite')")
}

# Step 2: Preserve the original Quartz PDF device on macOS.
# Cairo supplies editable PDF output on other platforms; font metrics may differ.
vector_pdf <- function(filename, width, height, family='Arial', bg='white', ...) {
  if (Sys.info()[['sysname']] == 'Darwin' && capabilities('aqua')) {
    grDevices::quartz(file=filename, type='pdf', width=width, height=height,
                     family=family, bg=bg, antialias=TRUE, ...)
  } else if (capabilities('cairo')) {
    grDevices::cairo_pdf(filename=filename, width=width, height=height,
                        family=family, bg=bg, ...)
  } else {
    stop('PDF export requires the Quartz or Cairo graphics device.')
  }
}

# Step 3: Export editable SVG, using fractional text sizing when supported.
vector_svg <- function(filename, width, height, family='Arial', bg='white', ...) {
  args <- list(filename=filename, width=width, height=height, bg=bg,
               system_fonts=list(sans=family))
  if ('fix_text_size' %in% names(formals(svglite::svglite))) {
    args$fix_text_size <- TRUE
  }
  do.call(svglite::svglite, c(args, list(...)))
}
