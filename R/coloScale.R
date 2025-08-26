#' Interpolate colors with values or classify values by colors
#'
#' This function is used to obtain a color scale from a set of values
#'
#' @param x A numeric vector
#' @param cols A vector of colors
#' @param center Whether to center the middle color value at 0
#' @param rescale01 Whether to keep values between 0 and 1 to a scale going from 0 to 1. Defaults to \code{FALSE}, meaning color values will go from min to max.
#' @param breaks A numeric vector. If breaks are given, value will be classified into each color provided the number of breaks equals the number of colors given. In other words, \code{length(breaks) == length(cols) + 1}
#'
#' @return A vector of interpolated colors
#'
#' @details
#' This function is a simple wrapper for the \code{\link{colorRamp}} function.
#'
#'
#' @export
#'
#'
coloScale<-function(
  x,
  cols=c("white","yellow","tomato3","darkred"),
  center=FALSE,
  rescale01=FALSE,
  breaks=NULL
){
  w<-which(is.na(x))
  if(any(w)){
    y<-x[-w]
  }else{
    y<-x
  }

  re<-function(a){
    if(any(w)){
      ans<-rep(NA,length(x))
      ans[which(!is.na(x))]<-a
      ans
    }else{
      a
    }
  }

  if(!is.null(breaks)){
    stopifnot((length(cols)+1)==length(breaks))
    return(re(cols[as.numeric(cut(y,breaks=breaks))]))
  }

  if(length(y)==1){
    colop<-colorRampPalette(cols)
    return(re(colop(y)))
  }
  if(class(y)=="character"){
    colop<-colorRampPalette(cols)
    color<-colop(length(unique(y)))
    return(re(color[match(y,unique(y))]))
  }else{
    if(all(y>=0 & y<=1) && rescale01){
      color<-rgb(colorRamp(cols)(y),maxColorValue=256)
      return(re(color))
    }else{
      if(center){
        m<-which.max(c(abs(min(y)),max(y)))
        sca<-0.5/ifelse(m==1,abs(min(y)),max(y))
        xx<-sca*y+0.5
        color<-rgb(colorRamp(cols)(xx),maxColorValue=256)
        return(re(color))
      }else{
        color<-rgb(colorRamp(cols)((y-min(y))/(max(y)-min(y))),maxColorValue=256)
        return(re(color))
      }
    }
  }

}
