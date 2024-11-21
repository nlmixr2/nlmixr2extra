if (!dir.exists("data")) {
  dir.create("data")
}

.in <- suppressWarnings(readLines("src/Makevars.in"))
if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win") ||
    (R.version$os == "linux-musl")) {
  file.out <- file("src/Makevars.win", "wb")
  writeLines(.in, file.out)
  close(file.out)
} else {
  file.out <- file("src/Makevars", "wb")
  writeLines(.in, file.out)
  close(file.out)
}
