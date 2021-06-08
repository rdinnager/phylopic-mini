library(imager)

unedited <- list.files("data/phylopic/mini_unedited", full.names = TRUE)

i <- 1L
for(i in seq_along(unedited)) {
  img <- imager::load.image(unedited[i])
  new_img <- cimg(array(0, dim = c(28, 28, 1, 4)))
  new_img[ , , 1, 4] <- img[ , , 1, 1] 
  
  imager::save.image(new_img, gsub("mini_unedited", "mini_transparent_bg", unedited[i]))
}