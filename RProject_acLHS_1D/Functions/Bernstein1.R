Fn.inv.Bernstein <- function(u, muestra){
  #
  # Entrada:
  #                 u = probabilidad para calcular cuantil
  #           muestra = vector de obervaciones univariadas x1,...,xn
  #
  x <- sort(muestra)
  n <- length(x)
  xm <- rep(0, n + 1)
  xm[1] <- x[1]
  xm[n + 1] <- x[n]
  xm[2:n] <- (x[1:(n-1)] + x[2:n])/2
  return(sum(xm * dbinom(0:n, n, u)))
}


Fn.Bernstein.aux <- function(u, xdada, muestra) Fn.inv.Bernstein(u, muestra) - xdada
Fn.Bernstein <- function(x, muestra) uniroot(Fn.Bernstein.aux, interval = c(0, 1),
                                             xdada = x, muestra = muestra, tol = 0.00001)$root


dFn.inv.Bernstein <- function(u, muestra){
  #
  # Entrada:
  #                 u = probabilidad para calcular cuantil
  #           muestra = vector de obervaciones univariadas x1,...,xn
  #
  x <- sort(muestra)
  n <- length(x)
  xm <- rep(0, n + 1)
  xm[1] <- x[1]
  xm[n + 1] <- x[n]
  xm[2:n] <- (x[1:(n-1)] + x[2:n])/2
  resultado <- x[n] * (u^(n - 1)) - x[1] * ((1 - u)^(n - 1))
  resultado <- resultado + sum(xm[2:n] * (dbinom(0:(n - 2), n - 1, u) - dbinom(1:(n - 1), n - 1, u)))
  resultado <- n * resultado
  return(resultado)
}
fn.Bernstein <- function(x, muestra) 1/dFn.inv.Bernstein(Fn.Bernstein(x, muestra), muestra)