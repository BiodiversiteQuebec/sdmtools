#' Infer habitat associations from a species sdm
#'
#' This function quantifies habitat associations by contrasting what is available in an area to what is use by the species
#'
#' @param raster A categorical SpatRaster with habitats.
#' @param sdm A SpatRaster with numerical values.
#' @param n Number of points to sample from the sdm. It can take a while to sample from sdm with high resolutions.
#' @param col Logical. Use colours associated with the raster if any? Default `TRUE`. if `FALSE`, no colours used.
#' @return A list with a ggplot object and a `data.frame` with the habitat association values.
#'
#' @details
#' This function compares the frequency of the different categories of habitats in the area to the frequencies of the habitats in the distribution of the species. First, it randomly samples points proportionally to the values in the sdm and calculates and calculates the relative proportions of the different habitats where the points fall. These relative proportions are illustrated in the top row of the graph. Second, it calculates the ratio between these proportions and the proportions of the different habitats in the habitat raster. Values over 1 means the habitat is selected, values below 1 means the habitat is avoided. The graph is shown on a log scale so that ratios of x have the same visual impact has a ratio of 1/x. Habitat categories where no random points are found (with a ratio of 0) are given the smallest ratio that is not 0.
#'
#'
#' @examples
#' \dontrun{
#'
#' library(terra)
#' library(sf)
#' library(rstac)
#' library(exactextractr)
#' library(rstac)
#' library(ggplot2)
#' library(geodata)
#' library(rmapshaper)
#'
#' # Downloads polygons using package geodata
#' can <- gadm("CAN", level = 1, path = getwd()) |> st_as_sf()
#' labrador <- ms_explode(can[can$NAME_1 %in% c("Newfoundland and Labrador"), ])
#' labrador <- labrador[which.max(st_area(labrador)), ] # keep Labarador
#' qc <- can[can$NAME_1 %in% c("Québec"),] |>
#'   rbind(labrador) |>
#'   ms_simplify(0.01)
#'
#'
#'
#' # https://www.cec.org/north-american-environmental-atlas/land-cover-30m-2020/
#'
#' category <- c("Temperate or sub-polar needleleaf forest",
#'              "Sub-polar taiga needleleaf forest",
#'              "Tropical or sub-tropical broadleaf evergreen forest",
#'              "Tropical or sub-tropical broadleaf deciduous forest",
#'              "Temperate or sub-polar broadleaf deciduous forest",
#'             "Mixed Forest",
#'              "Tropical or sub-tropical shrubland",
#'              "Temperate or sub-polar shrubland",
#'              "Tropical or sub-tropical grassland",
#'              "Temperate or sub-polar grassland",
#'             "Sub-polar or polar shrubland-lichen-moss",
#'              "Sub-polar or polar grassland-lichen-moss",
#'              "Sub-polar or polar barren-lichen-moss",
#'              "Wetland",
#'             "Cropland",
#'              "Barren lands",
#'              "Urban",
#'              "Water",
#'              "Snow and Ice")
#'
#' dat <- data.frame(code = c(0, seq_along(category)), category = c("No data", category))
#'
#' cec <- rast("/vsicurl/https://object-arbutus.cloud.computecanada.ca/bq-io/io/CEC_land_cover/NA_NALCMS_landcover_2020_30m.tif")
#'
#' levels(cec) <- dat
#'
#' cec <- crop(cec, st_transform(qc, st_crs(cec)))
#' raster <- aggregate(cec, 10, fun = "modal") |>
#'  crop(st_transform(qc, st_crs(cec)), mask = TRUE) |>
#'  droplevels()
#'
#' io <- stac("https://acer.biodiversite-quebec.ca/stac/")
#'
#'
#' coll <- "oiseaux"
#' #coll <- "sdm_emv_finance"
#'
#' urls <- io |>
#'   stac_search(collections = coll) |>
#'  post_request() |>
#'  items_fetch() |>
#'  _$features |>
#'  sapply(function(i){i$assets[[1]]$href})
#'
#' sdm <- rast(paste0('/vsicurl/', grep("(?=.*1900-2024)(?=.*setophaga_cer)", urls, perl = TRUE, value = TRUE)))
#'
#' x <- habitatAssociation(raster = raster, sdm = sdm, n = 10000, col = TRUE)
#' x
#'
#' }
#'
#' @import terra
#' @import sf
#' @import ggplot2
#'
#' @export
#'
#'
habitatAssociation <- function(
    raster,
    sdm,
    n = 10000,
    col = TRUE
) {

  if(nlyr(raster) > 1){
    raster <- raster[[1]]
  }

  sp <- spatSample(sdm, size = n, as.points = TRUE, method = "weights") |>
    st_as_sf() |>
    st_transform(st_crs(raster))

  v <- names(levels(raster)[[1]])[2]

  e1 <- extract(raster, sp)
  e1 <- table(e1[[v]]) |>
    as.data.frame() |>
    setNames(c(v, "count"))
  e1$p <- e1$count / sum(e1$count)

  e2 <- freq(raster)[, -1] |>
    setNames(c(v, "count"))

  e2[[v]] <- factor(e2[[v]], levels = levels(e1[[v]]))
  e2$p <- e2$count / sum(e2$count)


  e <- data.frame(category = e2[[v]], ratio = e1$p / e2$p) |>
    setNames(c(v, "ratio"))
  e$perc <- e1$p

  if(col){
    if(is.null(coltab(raster)[[1]])){
      warning("No coltab associated with the raster")
      e$col <- "forestgreen"
    } else {
      cols <- data.frame(code = coltab(raster)[[1]][,1], col = rgb(coltab(raster)[[1]][, 2:5] / 255))
      cols <- merge(cols, cats(raster)[[1]], all.x = TRUE)
      e$col <- cols$col[match(as.character(e[[v]]), cols[[v]])]
      #e$col[grep("Snow", e[[v]])] <- "lightsteelblue1"
    }
  } else {
    e$col <- "forestgreen"
  }

  e$ratio <- ifelse(e$ratio == 0, min(e$ratio[e$ratio > 0], na.rm = TRUE), e$ratio)
  e <- e[rev(order(e$ratio)),]
  e[[v]] <- factor(e[[v]], levels = e[[v]])

  g <- ggplot(data = e, aes(x = .data[[v]], y = ratio)) +
    geom_col(fill = e$col, width = 0.9) +
    #geom_point(aes(y = 1.5 * max(p)), fill = e$col, col = "white", shape = 21, size = e$width * 20, stroke = 0.0005) +
    geom_point(aes(y = 1.5 * max(ratio)), col = e$col, shape = 16, size = sqrt(e$perc) * 15) +
    #geom_col(fill = e$col, width = (e$width / max(e$width)) * 0.9) +
    scale_y_continuous(trans = "log", breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2, 5), limits = range(e$ratio) * c(1, 2)) +
    labs(x = "Habitat", y = "Utilisation / Disponibilité\n(évitement < 1 < sélection)\n", title = varnames(sdm)) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  list(habitats = e, plot = g)

}
