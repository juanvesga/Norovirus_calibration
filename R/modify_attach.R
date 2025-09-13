modify_attach <- function(pkg, include.only) {
	pkg <- as.character(substitute(pkg))

	if (!startsWith(pkg, "package:"))
		pkg <- paste0("package:", pkg)

	old <- if (pkg %in% search()) ls(pkg, all.names = TRUE)

	if (length(old))
		detach(pkg, character.only = TRUE)

	library(
		.rmpkg(pkg),
		include.only = unique(c(include.only, old)),
		character.only = TRUE
	)
}
