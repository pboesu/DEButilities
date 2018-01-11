
#' Scaled maturity volume at birth
#'
#' @param E_Hb maturity at birth
#' @param kap kappa
#' @param g energy investment ratio
#' @param p_Am maximum specific assimilation flux
#' @param v energy conductance
#' @param L_m maximum structural length
#'
#' @author Philipp Boersch-Supan
#'
#' @return Scaled maturity volume at birth v_Hb scalar
#' @export
#'
get_v_Hb <-  function(E_Hb, kap, g, p_Am, v, L_m){
  unname(E_Hb / ((1 - kap) * g * p_Am / v * L_m^3))
}
