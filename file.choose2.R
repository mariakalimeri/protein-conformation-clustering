file.choose2 <- function(...) {
    pathname <- NULL;
    tryCatch({
        pathname <- file.choose();
    }, error = function(ex) {
    })
    pathname;
}