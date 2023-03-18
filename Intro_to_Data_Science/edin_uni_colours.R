## University of Edinburgh Approved Colour Palettes
## source: https://uoe.sharepoint.com/sites/Brand/SitePages/Colours.aspx
## Author: Tim Old - tim.old@doctors.org.uk

## all colours specified seperately as both hexidecimal and RGB formats 
## default call to rgb() assume value of 0-1, need to set maxColorValue = 255 to match formatting
## from University of Edinburgh site
## Similarly, original hex values don't convert into R, inc because R format starts with "#" - 
## original from University of Edinburgh site commented out for reference, value used in assignment
## derived from call to rgb()


## Core colours - university red and university blue

edin_uni_red_hex <- "#D50032"
edin_uni_red_rgb <- rgb(193, 0, 67, maxColorValue = 255)

edin_uni_blue_hex <-  "#041E42"
edin_uni_blue_rgb <- rgb(4, 30, 66, maxColorValue = 255)

## Bright colours - bright pink, bright green, blue, bright orange,
## bright red, purple, bright yellow and bright blue

edin_bright_pink_hex <- "#D0006F"
edin_bright_pink_rgb <- rgb(208, 0, 111, maxColorValue = 255)

edin_bright_green_hex <- "#61BF1A"
edin_bright_green_rgb <- rgb(41, 188, 41, maxColorValue = 255)

edin_blue_hex <- "#0099AB"
edin_blue_rgb <- rgb(0, 114, 136, maxColorValue = 255)

edin_bright_orange_hex <- "#C25E03"
edin_bright_orange_rgb <- rgb(205, 90, 19, maxColorValue = 255)

edin_bright_red_hex <- "#AD033B"
edin_bright_red_rgb <- rgb(172, 0, 64, maxColorValue = 255)

edin_purple_hex <- "#830065"
edin_purple_rgb <- rgb(131, 0, 101, maxColorValue = 255)

edin_bright_yellow_hex <- "#F9A800"
edin_bright_yellow_rgb <- rgb(244, 170, 0, maxColorValue = 255)

edin_bright_blue_hex <- "#29C2DE"
edin_bright_blue_rgb <- rgb(0, 196, 223, maxColorValue = 255)

## all bright colours
edin_bright_cols_hex <- c(edin_bright_pink_hex, 
                          edin_bright_green_hex, 
                          edin_blue_hex, 
                          edin_bright_orange_hex,
                          edin_bright_red_hex, 
                          edin_purple_hex, 
                          edin_bright_yellow_hex, 
                          edin_bright_blue_hex)

edin_bright_cols_rgb <- c(edin_bright_pink_rgb, 
                          edin_bright_green_rgb, 
                          edin_blue_rgb, 
                          edin_bright_orange_rgb,
                          edin_bright_red_rgb, 
                          edin_purple_rgb, 
                          edin_bright_yellow_rgb, 
                          edin_bright_blue_rgb)


## Muted colours - muted brown, brown, muted yellow, muter turquoise,
## light blue, muted blue, dark green, muted pink

edin_muted_brown_hex <- "#704F45"
edin_muted_brown_rgb <- rgb(109, 79, 71, maxColorValue = 255)

edin_brown_hex <- "#692E1F"
edin_brown_rgb <- rgb(106, 51, 40, maxColorValue = 255)

edin_muted_yellow_hex <- "#949108"
edin_muted_yellow_rgb <- rgb(156, 154, 0, maxColorValue = 255)

edin_muted_turquoise_hex <- "#46877F"
edin_muted_turquoise_rgb <- rgb(69, 126, 129, maxColorValue = 255)

edin_light_blue_hex <- "#C6DBE9"
edin_light_blue_rgb <- rgb(194, 211, 223, maxColorValue = 255)

edin_muted_blue_hex <- "#004F71"
edin_muted_blue_rgb <- rgb(0, 79, 113, maxColorValue = 255)

edin_dark_green_hex <- "#154734"
edin_dark_green_rgb <- rgb(21, 71, 52, maxColorValue = 255)

edin_muted_pink_hex <- "BA8285"
edin_muted_pink_rgb <- rgb(184, 133, 141, maxColorValue = 255)

## all muted colours
edin_muted_cols_hex <- c(edin_muted_brown_hex, 
                         edin_brown_hex, 
                         edin_muted_yellow_hex, 
                         edin_muted_turquoise_hex, 
                         edin_light_blue_hex, 
                         edin_muted_blue_hex, 
                         edin_dark_green_hex, 
                         edin_muted_pink_hex)

edin_muted_cols_rgb <- c(edin_muted_brown_rgb, 
                         edin_brown_rgb, 
                         edin_muted_yellow_rgb, 
                         edin_muted_turquoise_rgb, 
                         edin_light_blue_rgb, 
                         edin_muted_blue_rgb, 
                         edin_dark_green_rgb, 
                         edin_muted_pink_rgb)

## Digital colours - university red, bright pink, dark green, muted brown,
## university blue, purple, jade, spruce grey,
## bright blue(2 - different to previous bright blue), burgundy, muted blue
## some repetition of previous - new only below

edin_jade_hex <- "487A7B"
edin_jade_rgb <- rgb(72, 122, 123, maxColorValue = 255)

edin_spruce_grey_hex <- "333F48"
edin_spruce_grey_rgb <- rgb(51, 63, 72, maxColorValue = 255)

edin_bright_blue_2_hex <- "007288"
edin_bright_blue_2_rgb <- rgb(0, 114, 136, maxColorValue = 255)

edin_burgundy_hex <- "A50034"
edin_burgundy_rgb <- rgb(165, 0, 52, maxColorValue = 255)

## all digital colours

edin_dig_cols_hex <- c(edin_uni_red_hex,
                       edin_bright_pink_hex, 
                       edin_dark_green_hex, 
                       edin_muted_brown_hex, 
                       edin_uni_blue_hex,
                       edin_purple_hex, 
                       edin_jade_hex, 
                       edin_spruce_grey_hex, 
                       edin_bright_blue_2_hex, 
                       edin_burgundy_hex, 
                       edin_muted_blue_hex)

edin_dig_cols_rgb <- c(edin_uni_red_rgb,
                       edin_bright_pink_rgb, 
                       edin_dark_green_rgb, 
                       edin_muted_brown_rgb, 
                       edin_uni_blue_rgb,
                       edin_purple_rgb, 
                       edin_jade_rgb, 
                       edin_spruce_grey_rgb, 
                       edin_bright_blue_2_rgb, 
                       edin_burgundy_rgb, 
                       edin_muted_blue_rgb)