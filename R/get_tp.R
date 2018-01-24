
#' Gets scaled age at puberty
#'
#' @param p 5-vector with parameters: g, k, l_T, v_H^b, v_H^p
#' @param f optional scalar with functional response (default f = 1)
#' @param lb0 optional scalar with scaled length at birth or optional 2-vector with scaled length, l, and scaled maturity, vH, for a juvenile that is now exposed to f, but previously at another f
#'
#' @return [tp, tb, lp, lb, info]
#' @export
#'
#' @examples get_tp(c(.5, .1, .1, .01, .2))
get_tp <- function(p, f = 1, lb0 = NULL){
## get_tp
# Gets scaled age at puberty

##
  #function [tp, tb, lp, lb, info] = get_tp(p, f, lb0)
  # created at 2008/06/04 by Bas Kooijman,
  # modified 2014/03/04 Starrlight Augustine, 2015/01/18 Bas Kooijman

  ## Syntax
  # [tp, tb, lp, lb, info] = <../get_tp.m *get_tp*>(p, f, lb0)

  ## Description
  # Obtains scaled age at puberty.
  # Food density is assumed to be constant.
  # Multiply the result with the somatic maintenance rate coefficient to arrive at age at puberty.
  #
  # Input
  #
  # * p: 5-vector with parameters: g, k, l_T, v_H^b, v_H^p
  # * f: optional scalar with functional response (default f = 1)
  # * lb0: optional scalar with scaled length at birth
  #
  #      or optional 2-vector with scaled length, l, and scaled maturity, vH
  #      for a juvenile that is now exposed to f, but previously at another f
  #
  # Output
  #
  # * tp: scaled with age at puberty \tau_p = a_p k_M
  #
  #      if length(lb0)==2, tp is the scaled time till puberty
  #
  # * tb: scaled with age at birth \tau_b = a_b k_M
  # * lp: scaler length at puberty
  # * lb: scaler length at birth
  # * info: indicator equals 1 if successful

  ## Remarks
  #  Function get_tp_foetus does the same for foetal development; the result depends on embryonal development.

  ## Example of use
  # get_tp([.5, .1, .1, .01, .2])

  #  unpack pars
  g   = p[1]; # energy investment ratio
  k   = p[2]; # k_J/ k_M, ratio of maturity and somatic maintenance rate coeff
  lT  = p[3]; # scaled heating length {p_T}/[p_M]Lm
  vHb = p[4]; # v_H^b = U_H^b g^2 kM^3/ (1 - kap) v^2; U_H^b = M_H^b/ {J_EAm} = E_H^b/ {p_Am}
  vHp = p[5]; # v_H^p = U_H^p g^2 kM^3/ (1 - kap) v^2; U_H^p = M_H^p/ {J_EAm} = E_H^p/ {p_Am}

  #if (is.null(lb0)) lb0 = NA;

  if (k == 1 && f * (f - lT)^2 > vHp * k){
    lb = vHb^(1/3);
    tb = get_tb(p =c(1, 2, 4), f, lb);
    lp = vHp^(1/3);
    li = f - lT;
    rB = 1 / 3/ (1 + f/g);
    tp = tb + log((li - lb)/ (li - lp))/ rB;
    info = 1;
  } else {
    if (f * (f - lT)^2 <= vHp * k){ # reproduction is not possible
      pars_lb = c(1, 2, 4);

      tb_lb_info = get_tb (pars_lb, f, lb0);
      tb = tb_lb_info["tb"]
      lb = tb_lb_info["lb"]
      info = tb_lb_info["info"]
      tp = 1e20; # tau_p is never reached
      lp = 1;    # lp is nerver reached
    } else { # reproduction is possible
      li = f - lT;
      irB = 3 * (1 + f/ g); # k_M/ r_B
      lp_lb_info <- get_lp(p, f , lb0)
      lp = lp_lb_info["lp"]
      lb = lp_lb_info["lb"]
      info = lp_lb_info["info"]
      #[lp, lb, info] = get_lp(p, f, lb0);
        if (length(lb0) != 2){ # lb0 = l_b
        tb = get_tb(c(g, k, vHb), f, lb)["tb"];
        tp = tb + irB * log((li - lb)/ (li - lp));
        } else {# lb0 = l and t for a juvenile
        tb = NaN;
        l = lb0[1];
        tp = irB * log((li - l)/ (li - lp));
        }}}

  if (!is.double(tp) || tp < 0){ # tp must be real and positive
    info = 0;
  }

return(c(tp = unname(tp), tb = unname(tb), lp = unname(lp), lb = unname(lb), info = unname(info)))
}
