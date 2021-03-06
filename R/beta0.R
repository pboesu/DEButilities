#' particular incomplete beta function
#'
#' \deqn{B_x1(4/3,0) - B_x0(4/3,0) = \int_x0^x1 t^(4/3-1) (1-t)^(-1) dt}{B_x1(4/3,0) - B_x0(4/3,0) = int_x0^x1 t^(4/3-1) (1-t)^(-1) dt}
#'
#' @param x0 scalar with lower boundary for integration
#' @param x1 scalar with upper boundary for integration
#'
#' @return scalar with particular incomple beta function
#' @export
#'
#' @examples
#' beta0(0.1, 0.2)
#'
beta0 = function(x0,x1){
  #  created 2000/08/16 by Bas Kooijman; modified 2011/04/10

  ## Syntax
  # f = <../beta.m *beta0*> (x0,x1)

  ## Description
  #  particular incomplete beta function:
    #   B_x1(4/3,0) - B_x0(4/3,0) = \int_x0^x1 t^(4/3-1) (1-t)^(-1) dt

  # Input
  #
  # * x0: scalar with lower boundary for integration
  # * x1: scalar with uper boundary for integration
  #
  # Output
  #
  # * f: scalar with particular incomple beta function

  ## Remarks
  # See also <../lib/misc/beta_34_0 *beta_34_0*>

    ## Example of use
  # beta0(0.1, 0.2)

  if (x0 < 0 || x0 >= 1 || x1 < 0 || x1 >= 1){
    print('Warning from beta0: argument values outside (0,1) \n')
    f = NA
  }

  n0 = length(x0)
  n1 = length(x1)
  if (n0 != n1 && n0 != 1 && n1 != 1) {
    print('Warning from beta0: argument sizes do not match \n')
    f = NA
  }

  x03 = x0 ^ (1/ 3)
  x13 = x1 ^ (1/ 3)
  a3 = sqrt(3)
  f1 = - 3 * x13 + a3 * atan((1 + 2 * x13)/ a3) - log(as.complex(x13 - 1)) + log(as.complex(1 + x13 + x13 ^ 2))/ 2
  f0 = - 3 * x03 + a3 * atan((1 + 2 * x03)/ a3) - log(as.complex(x03 - 1)) + log(as.complex(1 + x03 + x03 ^ 2))/ 2
  f = Re(f1 - f0)

  return(f)
}
