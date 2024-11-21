if (!dir.exists("data")) {
  dir.create("data")
}

.in <- suppressWarnings(readLines("src/Makevars.in"))
if (.Platform$OS.type == "windows") {
  .makevars <- file("src/Makevars.win", "wb")
  .i <- "I"
} else {
  .makevars <- file("src/Makevars", "wb")
  if (file.exists("/etc/os-release")) {
    .os <- readLines("/etc/os-release")
    if (any(grepl("Pop!_OS", .os, fixed=TRUE))) {
      .i <- "isystem"
    } else {
      .i <- "I"
    }
  } else {
    .i <- "I"
  }
}

writeLines(gsub("@ISYSTEM@", .i, .in),
           .makevars)
close(.makevars)
