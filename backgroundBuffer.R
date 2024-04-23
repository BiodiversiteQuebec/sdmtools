#' Add background points outside a buffer around occurrences
#'
#' This function generates a buffer around occurrence points and samples points outside the buffer to create a background where the species will be interpreted as being absent outside this buffer.
#' 
#' @param occs An sf object with occurrence points.
#' @param region An sf object with the study region.
#' @param n The number of background points to sample.
#' @param dist The distance for buffering around occurrence points. The default is 500 km. The value has to be given in the units of the crs (i.e. often in meters).
#' @param convex Logical. If `TRUE`, the buffer is computed on a convex hull around occurrences points. Default is `FALSE`. A concave hull could also be used with the package concaveman (but not used at the moment).  
#' @param prop Logical. If `TRUE`, the number of points sampled will be such that the density of points in the area outside of the buffer will be equal to the density of occurrence points in the convex hull around occurrences or within the buffer area. Default `FALSE`. 
#' @param dens_buffer Logical. If `TRUE`, the density of background points will be equal to the density of occurrence points inside the buffer area. If `FALSE`, the density will be equal to the density of occurrences points within the convex hull around occurrences.
#' @param nmax Maximum number of points to sample when `prop = FALSE`
#' 
#' @return An sf object with the background points. Returns `NULL` if there is no area outside of the buffer used for the study area.
#' 
#' @details
#' This function is used to add an artifical effort where we are relatively certain the target species is not found. By adding background points outside of the chosen buffer, it can help in reducing the lkelihood that the sdm model predicts high presence/abundance of the species outside of its actual range. The use of this function relies on the assumption that the species is not found outside of the buffer area. The function is meant to be used in conjunction with a target-group approach where background points are also generated using a taxonomic group with similar bias. Both types of background points can be added. The target-group background points ensure that the sampling bias is taken into account within the buffer area and the background points generated through the function here ensure that predictions do not extend much outside of the buffer area.  
#' 
#' 
#' @examples
#' library(sf)
#' 
#' @import sf
#' 
#' @export
#' 
backgroundBuffer <- function(
    occs, 
    region, 
    n = 1000, 
    dist = 500000, 
    convex = FALSE,
    prop = FALSE, 
    dens_buffer = FALSE,
    nmax = NULL
) {
  x <- st_union(occs) |> st_as_sf()   
  ### When the number of occs becomes larger, using concaveman could be much faster than st_buffer
  if(convex) {
    buffer <- st_buffer(st_convex_hull(x), dist = dist)
  } else {
    buffer <- st_buffer(st_cast(x, "POINT"), dist = dist) |> st_union() |> st_as_sf()
  }
  outside <- st_difference(st_geometry(region), buffer)
  if(is.null(length(outside)) || length(outside) == 0){
    return(NULL)
  }
  if(prop) {
    np <- nrow(st_cast(x, "POINT")) # if occs is unioned, we explode the points to get the number of points
    if(dens_buffer) {
      a <- as.numeric(st_area(st_intersection(st_geometry(region), buffer)))
    } else {
      a <- as.numeric(st_area(st_convex_hull(x)))
    }
    ns <- round((np / a) * as.numeric(st_area(outside)))
    if(ns < 1) { 
      ns <- 1 
    }
    if(!is.null(nmax)){
      ns <- min(c(ns, nmax))
    }
    s <- st_sample(outside, size = ns, type = "random")
  } else {
    s <- st_sample(outside, size = n, type = "random")
  }
  st_as_sf(s) 
}



library(sf)
library(geodata)
library(rmapshaper)
library(rgbif)

epsg <- 6622

can <- gadm('CAN', level = 1, path = getwd()) |> st_as_sf()
qc <- can[can$NAME_1 == "Québec", ]
qc <- ms_simplify(qc, 0.001) 
qc <- st_transform(qc, crs = epsg)
region <- qc

xx <- occ_data(scientificName = "Fagus grandifolia", recordedBy = NULL, hasCoordinate = TRUE, limit = 5500, stateProvince = "Québec")
x <- st_as_sf(xx$data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
x <- st_transform(x, crs = epsg)

occs <- x[region, ]

plot(st_geometry(region))
plot(st_geometry(occs), add = TRUE)

par(mfrow=c(2,3), mar = c(0, 0, 0, 0))

# fixed number of background points
plot(st_geometry(region))
plot(st_geometry(occs), add = TRUE)
background <- backgroundBuffer(occs, region, n = 500) 
plot(background, col = "red", add = TRUE)

# same as the previous one with a 100 km buffer
plot(st_geometry(region))
plot(st_geometry(occs), add = TRUE)
background <- backgroundBuffer(occs, region, n = 500, dist = 100000) 
plot(background, col = "red", add = TRUE)

# same as the previous one, but buffer around the convex hull
plot(st_geometry(region))
plot(st_geometry(occs), add = TRUE)
background <- backgroundBuffer(occs, region, n = 500, dist = 100000, convex = TRUE) 
plot(background, col = "red", add = TRUE)

# density of background points equals to the density within the convex hull around occurrences
plot(st_geometry(region))
plot(st_geometry(occs), add = TRUE)
background <- backgroundBuffer(occs, region, prop = TRUE) 
plot(background, col = "red", add = TRUE)

# same as the previous one, but with a maximum of 200 points used
plot(st_geometry(region))
plot(st_geometry(occs), add = TRUE)
background <- backgroundBuffer(occs, region, prop = TRUE, nmax = 200) 
plot(background, col = "red", add = TRUE)

# density of background points equals to the density inside the buffer area
plot(st_geometry(region))
plot(st_geometry(occs), add = TRUE)
background <- backgroundBuffer(occs, region, prop = TRUE, dens_buffer = TRUE) 
plot(background, col = "red", add = TRUE)










