# Define global variables to satisfy R CMD check constraints regarding ggplot2 aesthetics
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("X", "Y", "Expression", "spatial_effect_size", "adjusted_p_value", "Significance"))
}
