############################################################################################
## Content: R script to calculate the drift of a moving object subject to an external force.
## Date: 14/11/2017
## Author: Karel De Vogeleer
############################################################################################

# Degree to radian convertors
deg2rad <- function(deg) {deg/180*pi}
rad2deg <- function(rad) {rad*180/pi}

# Converts an arbitrary angle to an equivalent angle: 0 <= angle < 360 
# Warning: this is a recursive angle, cannot scale infinitely
normalize_angle <- function(hdg) {
  if(hdg >= 0 & hdg < 360)
    return(hdg)
  if(hdg < 0)
    return(normalize_angle(hdg+360))
  if(hdg >= 360)
    return(normalize_angle(hdg-360))
}

# Calculates the drift of an object based on the cosine rule
calculate_drift_cosine <- function(object, force) {
  force_original <- force
  object_original <- object
  
  output <- c()
  force$heading <- normalize_angle(force$heading-object$heading)
  angle <- abs(180-force$heading)
    
  B <- -2*cos(deg2rad(angle))*force$velocity
  C <- force$velocity^2-object$velocity^2
  output$indicated_velocity <- Re(polyroot(c(C,B,1))[1])
  output$correction <- rad2deg(acos((force$velocity^2-(output$indicated_velocity^2+object$velocity^2))/(-2*output$indicated_velocity*object$velocity)))
  
  if(force$heading > 180)
    output$correction <- output$correction * -1

  output$Tc <- object_original$heading
  output$Tt <- normalize_angle(object_original$heading + output$correction)
    
  return(output)
} 

# Calculates the drift of an object based on the theorem of Pythogoras
calculate_drift_pythagoras <- function(object, force) {
  force_original <- force
  object_original <- object
  
  output <- c()
  force$heading <- normalize_angle(force$heading-object$heading)
  angle <- abs(180-force$heading)
  
  output$correction <- rad2deg(asin(sin(deg2rad(angle))*force$velocity/object$velocity))
  output$indicated_velocity <- cos(deg2rad(angle))*force$velocity+cos(deg2rad(output$correction))*object$velocity
  
  if(force$heading > 180)
    output$correction <- output$correction * -1
  
  output$Tc <- object_original$heading
  output$Tt <- normalize_angle(object_original$heading + output$correction)
  
  return(output)
}

# Units:
# "heading" is assumed to be in degrees
# "velocity" can be anything (e.g., m/s or kts) as long as the same unit is used consistently

# Define the object
object <- c()
object$heading <- 27
object$velocity <- 95

# Define the force
#  The force heading is the source of the force, not the force direction. This is common for wind calculations, for example.  
force <- c()
force$heading <- 243
force$velocity <- 15

# Calulate drift via cosine and Pythagoras
# Both outcomes should be the same
print(calculate_drift_pythagoras(object,force))
print(calculate_drift_cosine(object,force))
