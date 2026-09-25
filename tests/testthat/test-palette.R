test_that("magreport_palette returns distinct colours", {
  expect_length(magreport_palette(0), 0)
  expect_identical(magreport_palette(3), magreport_palette(12)[1:3])
  expect_equal(anyDuplicated(magreport_palette(12)), 0)
  expect_equal(anyDuplicated(magreport_palette(20)), 0)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", magreport_palette(20))))
})

test_that("quality_colours covers every quality level", {
  expect_named(quality_colours(), quality_levels())
  expect_true(all(quality_colours() %in% magreport_colours()))
})

test_that("scales and theme work with ggplot2", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(x = 1:3, y = 1:3, g = c("a", "b", "c"))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = g, fill = y)) +
    ggplot2::geom_point() +
    scale_colour_magreport() +
    scale_fill_magreport_c() +
    theme_magreport()

  built <- ggplot2::ggplot_build(p)
  expect_setequal(unique(built$data[[1]]$colour), magreport_palette(3))
  expect_s3_class(theme_magreport(), "theme")
  expect_s3_class(scale_fill_magreport(), "ScaleDiscrete")
})
