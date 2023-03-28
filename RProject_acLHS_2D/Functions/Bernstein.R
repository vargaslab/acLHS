#----------------------------
# Cópula Bernstein bivariada
#----------------------------

a_leer_primero <- function() {
    cat(noquote("Los datos bivariados (x1,y1),...,(xn,yn) deben estar capturados en una matriz"), "\n")
    cat(noquote("de n x 2 con el nombre <muestra>"), "\n")
    cat(noquote("Genere la matriz de la c'opula emp'irica mediante la siguiente instrucci'on:"), "\n")
    cat(noquote("matriz.copem <- genmat.copem(muestra)"), "\n")
    cat(noquote("Hecho lo anterior puede proceder a los c'alculos."), "\n")
    cat(noquote("Fin instructivo."), "\n")
}

Bpol.fun <- function(x, funcion, orden) sum(funcion((0:orden)/orden) * dbinom(0:orden, orden, x))

Bpol.fun.emp <- function(x, valores.emp) sum(valores.emp * dbinom(0:(length(valores.emp) - 1), length(valores.emp) - 1, x))

Bpol.valores <- function(valores, funcion, orden) {
    vec.valores <- rep(0, length(valores))
    for (j in 1:(length(valores))) {
        vec.valores[j] <- Bpol.fun(valores[j], funcion, orden)
    }
    return(vec.valores)
}

Bpol.valores.emp <- function(valores, valores.emp) {
    vec.valores <- rep(0, length(valores))
    for (j in 1:(length(valores))) {
        vec.valores[j] <- Bpol.fun.emp(valores[j], valores.emp)
    }
    return(vec.valores)
}

copula.Bernshtein.emp <- function(u, v) sum(matriz.copem * (dbinom(0:(dim(matriz.copem)[1] - 1), dim(matriz.copem)[1] - 1, u) %*% t(dbinom(0:(dim(matriz.copem)[1] - 1), dim(matriz.copem)[1] - 1, v))))

cv.du <- function(u, v) (dim(matriz.copem)[1] - 1) * sum(matriz.copem * ((dbinom(-1:(dim(matriz.copem)[1] - 2), dim(matriz.copem)[1] - 2, u) - dbinom(0:(dim(matriz.copem)[1] - 1), 
                         dim(matriz.copem)[1] - 2, u) * c(-1, rep(1, dim(matriz.copem)[1] - 1))) %*% t(dbinom(0:(dim(matriz.copem)[1] - 1), dim(matriz.copem)[1] - 1, v))))

cv.du.aux <- function(v, ua.vec) cv.du(ua.vec[1], v) - ua.vec[2]

cv.du.inv <- function(u, a) uniroot(cv.du.aux, interval = c(0, 1), ua.vec = c(u, a), tol = tolerancia)$root

dcopula.Bernshtein.emp <- function(u, v) ((dim(matriz.copem)[1] - 1)^2) * sum(matriz.copem * ((dbinom(-1:(dim(matriz.copem)[1] - 2), dim(matriz.copem)[1] - 2, u) - dbinom(0:(dim(matriz.copem)[1] - 
                                          1), dim(matriz.copem)[1] - 2, u) * c(-1, rep(1, dim(matriz.copem)[1] - 1))) %*% t(dbinom(-1:(dim(matriz.copem)[1] - 2), dim(matriz.copem)[1] - 
                                          2, v) - dbinom(0:(dim(matriz.copem)[1] - 1), dim(matriz.copem)[1] - 2, v) * c(-1, rep(1, dim(matriz.copem)[1] - 1)))))

densidad.Bernshtein.emp <- function(x, y) dcopula.Bernshtein.emp(Fn.Bernshtein(x), Gn.Bernshtein(y))*fn.Bernshtein(x)*gn.Bernshtein(y)

densidad.Bernshtein.emp.ydadox <- function(y, x) dcopula.Bernshtein.emp(Fn.Bernshtein(x), Gn.Bernshtein(y))*gn.Bernshtein(y)

dep.schweizer <- function (mat.xy){
    copem <- genmat.copem(mat.xy)
    n <- dim(mat.xy)[1]
    mat1 <- matrix(seq(from = 0, to = 1, length = (n + 1)), nrow = (n + 
        1), ncol = 1)
    mat2 <- matrix(seq(from = 0, to = 1, length = (n + 1)), ncol = (n + 
        1), nrow = 1)
    copula.pi <- mat1 %*% mat2
    sigma.schweizer <- (12/(n^2)) * sum(abs(copem - copula.pi))
    return(sigma.schweizer)
}

depurar.xy <- function(mat.xy){
    n <- dim(mat.xy)[1]
    matriz.ord <- mat.xy
    orden <- order(mat.xy[, 1])
    for (i in 1:n) {
        matriz.ord[i, ] <- mat.xy[orden[i], ]
    }
    i <- 1
    rango.bloque <- which(matriz.ord[, 1] == matriz.ord[i, 1])
    matriz.dep <- matrix(c(matriz.ord[i, 1], median(matriz.ord[rango.bloque, 2])), ncol = 2, nrow = 1)
    i <- max(rango.bloque) + 1
    while (i <= n) {
        rango.bloque <- which(matriz.ord[, 1] == matriz.ord[i, 1])
        matriz.dep <- rbind(matriz.dep, c(matriz.ord[i, 1], median(matriz.ord[rango.bloque, 2])))
        i <- max(rango.bloque) + 1
    }
    mat.xy <- matriz.dep[, 2:1]
    n <- dim(mat.xy)[1]
    matriz.ord <- mat.xy
    orden <- order(mat.xy[, 1])
    for (i in 1:n) {
        matriz.ord[i, ] <- mat.xy[orden[i], ]
    }
    i <- 1
    rango.bloque <- which(matriz.ord[, 1] == matriz.ord[i, 1])
    matriz.dep <- matrix(c(matriz.ord[i, 1], median(matriz.ord[rango.bloque, 2])), ncol = 2, nrow = 1)
    i <- max(rango.bloque) + 1
    while (i <= n) {
        rango.bloque <- which(matriz.ord[, 1] == matriz.ord[i, 1])
        matriz.dep <- rbind(matriz.dep, c(matriz.ord[i, 1], median(matriz.ord[rango.bloque, 2])))
        i <- max(rango.bloque) + 1
    }
    mat.xy <- matriz.dep[, 2:1]
    n <- dim(mat.xy)[1]
    matriz.ord <- mat.xy
    orden <- order(mat.xy[, 1])
    for (i in 1:n) {
        matriz.ord[i, ] <- mat.xy[orden[i], ]
    }
    return(matriz.ord)
}

dFn.inv.Bernshtein <- function(u, valores.emp){
    x <- sort(valores.emp)
    n <- length(x)
    xm <- rep(0, n + 1)
    for (j in 2:n) {
        xm[j] <- (x[j - 1] + x[j])/2
    }
    xm[1] <- x[1]
    xm[n + 1] <- x[n]
    resultado <- x[n] * (u^(n - 1)) - x[1] * ((1 - u)^(n - 1))
    resultado <- resultado + sum(xm[2:n] * (dbinom(0:(n - 2), n - 1, u) - dbinom(1:(n - 1), n - 1, u)))
    resultado <- n * resultado
    return(resultado)
}

diag.copem <- function(matcopem){
  m <- ncol(matcopem) - 1
  return(cbind((0:m)/m, diag(matcopem)))
}

estandarizar <- function(muestra) apply(muestra, 2, rank)/nrow(muestra) # calcular muestra estandarizada

fn.Bernshtein <- function(x) 1/dFn.inv.Bernshtein(Fn.Bernshtein(x), muestra[, 1])

Fn.Bernshtein <- function(x) uniroot(Fn.Bernshtein.aux, interval = c(0, 1), xdada = x, tol = tolerancia, extendInt="yes", trace=2)$root

Fn.Bernshtein.aux <- function(u, xdada) Fn.inv.Bernshtein(u, muestra[, 1]) - xdada

Fn.emp <- function(x, datos) mean(datos <= x)

Fn.inv.Bernshtein <- function(u, valores.emp){
    x <- sort(valores.emp)
    n <- length(x)
    xm <- rep(0, n + 1)
    for (j in 2:n) {
        xm[j] <- (x[j - 1] + x[j])/2
    }
    xm[1] <- x[1]
    xm[n + 1] <- x[n]
    return(sum(xm * dbinom(0:n, n, u)))
}

genmat.copem <- function(mat.xy){
    n <- dim(mat.xy)[1]
    mat.copem <- matrix(0, ncol = (n + 1), nrow = (n + 1))
    mat.xyord <- mat.xy
    orden <- order(mat.xy[, 1])
    for (i in 1:n) {
        mat.xyord[i, ] <- mat.xy[orden[i], ]
    }
    mat.copem[n + 1, ] <- (0:n)/n
    y.ord <- sort(mat.xyord[, 2])
    for (i in 1:(n - 1)) {
        columna <- (((mat.xyord[, 2][i] <= y.ord)) * 1)/n
        mat.copem[i + 1, ] <- mat.copem[i, ] + c(0, columna)
    }
    return(mat.copem)
}

genmat.copem.Bernshtein <- function(u.vec, v.vec){
    copula.B <- matrix(0, nrow = length(u.vec), ncol = length(v.vec))
    for (i in 1:(length(u.vec))) {
        for (j in 1:(length(v.vec))) {
            copula.B[i, j] <- copula.Bernshtein.emp(u.vec[i], v.vec[j])
        }
    }
    return(list(u = u.vec, v = v.vec, copemB = copula.B))
}

genmat.dcopem.Bernshtein <- function(u.vec, v.vec){
    dcopula.B <- matrix(0, nrow = length(u.vec), ncol = length(v.vec))
    for (i in 1:(length(u.vec))) {
        for (j in 1:(length(v.vec))) {
            dcopula.B[i, j] <- dcopula.Bernshtein.emp(u.vec[i], 
                v.vec[j])
        }
    }
    return(list(u = u.vec, v = v.vec, dcopemB = dcopula.B))
}

genmat.densidad.Bernshtein <- function(x.vec, y.vec){
    densidad <- matrix(0, nrow = length(x.vec), ncol = length(y.vec))
    for (i in 1:(length(x.vec))) {
        for (j in 1:(length(y.vec))) {
            densidad[i, j] <- densidad.Bernshtein.emp(x.vec[i], 
                y.vec[j])
        }
    }
    return(list(x = x.vec, y = y.vec, densidad = densidad))
}

genmat.diagem <- function(matriz){
#
# Input: matriz de (n x 2) cuyos renglones son observaciones del vector aleatorio (X,Y).
# Output: matriz de (n+1 x 2) con los valores (u,dn(u)) de la diagonal empirica
#
    n <- dim(matriz)[1]     #  tamanio de muestra
    tau <- rep(0, n-1)      #  declarar subvector de trayectoria de tau(1) a tau(n-1)
    x <- matriz[, 1]
    y <- matriz[, 2]
    x.orden <- rank(x)      #  vector de orden en los valores de x
    y.ordx <- y             #  declarando vector de y ordenado en x
    for(j in 1:n){          #  ordenando y de acuerdo a x
        y.ordx[x.orden[j]] <- y[j]
    }
    y <- sort(y)            #  ordenando y de menor a mayor
    for(j in 1:(n-1)){      #  calculando vector tau[1:n-1]
        acumulador <- 0
        for(k in 1:j){
            acumulador <- acumulador + 1*(y.ordx[k] <= y[j])
        }
        tau[j] <- acumulador
    }
    tau <- c(0, tau, n)/n   # vector completo tau[0:n]
    mat.diag <- matrix(0, ncol = 2, nrow = (n+1))
    mat.diag[ , 1] <- (0:n)/n
    mat.diag[ , 2] <- tau
    return(mat.diag)
}

gn.Bernshtein <- function(x) 1/dFn.inv.Bernshtein(Gn.Bernshtein(x), muestra[, 2])

Gn.Bernshtein <- function(x) uniroot(Gn.Bernshtein.aux, interval = c(0, 1), xdada = x, tol = tolerancia, extendInt="yes", trace=2)$root

Gn.Bernshtein.aux <- function(u, xdada) Fn.inv.Bernshtein(u, muestra[, 2]) - xdada

grafica.cotas.diagonal <- function(titulo){
#
# Genera una plantilla para graficar u versus la
# secci?n diagonal d(u), agregando cotas de FH y
# la diagonal de Pi
#
# Input: t?tulo del gr?fico (entrecomillado)
#
  plot(c(0, 1),c(0, 1),type = "n", main = titulo, xlab = "u", ylab = "diag(u)")
  lines(c(0, 1),c(0, 1), col = "green")
  lines(c(0, 0.5), c(0, 0), col="green")
  lines(c(0.5, 1), c(0, 1), col = "green")
  u <- seq(from = 0, to = 1, length = 1000)
  lines(u, u^2, col = "orange")
}

regresion <- function(x, cuantil) Fn.inv.Bernshtein(regresion.copulaB(Fn.Bernshtein(x), cuantil), muestra[, 2])

regresion.copulaB <- function(u, cuantil) uniroot(cv.du.aux, interval = c(0, 1), ua.vec = c(u, cuantil), tol = tolerancia)$root

simula.Bernshtein <- function(tam.muestra){
    sim.copula <- simula.copula.Bernshtein(tam.muestra)
    x <- sapply(sim.copula[, 1], Fn.inv.Bernshtein, valores.emp = muestra[, 1])
    y <- sapply(sim.copula[, 2], Fn.inv.Bernshtein, valores.emp = muestra[, 2])
    simulaciones <- cbind(x, y)
    return(list(sim.xy = simulaciones, sim.copula = sim.copula))
}

simula.Bernshtein.condicional <- function(x, tam.muestra){
    uu <- runif(tam.muestra)
    Fx <- Fn.Bernshtein(x)
    v <- sapply(uu, cv.du.inv, u = Fx)
    y <- sapply(v, Fn.inv.Bernshtein, valores.emp = muestra[, 2])
    return(y)
}

simula.copula.Bernshtein <- function(tam.muestra){
    uu <- runif(tam.muestra)
    tt <- runif(tam.muestra)
    vv <- mapply(cv.du.inv, u = uu, a = tt)
    return(cbind(uu, vv))
}

simula.copula.unif3 <- function(tam.muestra, theta){
    n <- tam.muestra
    m <- matrix(0, ncol = 2, nrow = n)
    m[, 1] <- runif(n)
    a <- rbinom(n, 1, theta)
    b <- 2 * rbinom(n, 1, 1/2) - 1
    m[, 2] <- (m[, 1] <= (1/3)) * (runif(n, 0, 2/3) * a + runif(n, 2/3, 1) * (1 - a)) + (((1/3) < m[, 1]) & (m[, 1] <= (2/3))) * (runif(n, 1/3, 1) * a + runif(n, 0, 1/3) * (1 - a)) + ((2/3) < m[, 1]) * (runif(n, 1/3, 2/3) + a * b/3)
    return(m)
}

transf.muestra.copula <- function (){
    u <- sapply(muestra[, 1], Fn.Bernshtein)
    v <- sapply(muestra[, 2], Gn.Bernshtein)
    return(cbind(u, v))
}

tolerancia <- 0.00001


