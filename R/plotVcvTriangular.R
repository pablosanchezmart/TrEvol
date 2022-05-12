plotVcvTriangular <- function (corr, corrProp, col = NULL, col.lim = NULL, bg = "white", title = "", add = FALSE, diag = TRUE, outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL,  tl.pos = "ld",
                               tl.cex = 1, tl.col = "black", tl.offset = 0.4, tl.srt = 90, cl.pos = "b", cl.length = NULL, cl.cex = 0.8, cl.ratio = 0.15, cl.align.text = "c", cl.offset = 0.5, type = "lower")
{


  ### FUNCTIONS ####

  draw_grid = function(coords, fg) {
    symbols(coords, add = TRUE, inches = FALSE, fg = fg, bg = NA,
            rectangles = matrix(1, nrow = nrow(coords), ncol = 2))
  }

  apply_mat_filter = function(mat) {
    x = matrix(1:n * m, nrow = n, ncol = m)
    switch(type, upper = mat[row(x) > col(x)] <- Inf, lower = mat[row(x) <
                                                                    col(x)] <- Inf)
    if (!diag) {
      diag(mat) = Inf
    }
    return(mat)
  }

  getPos.Dat = function(mat) {
    tmp = apply_mat_filter(mat)
    Dat = tmp[is.finite(tmp)]
    ind = which(is.finite(tmp), arr.ind = TRUE)
    Pos = ind
    Pos[, 1] = ind[, 2]
    Pos[, 2] = -ind[, 1] + 1 + n
    PosName = ind
    PosName[, 1] = colnames(mat)[ind[, 2]]
    PosName[, 2] = rownames(mat)[ind[, 1]]
    return(list(Pos, Dat, PosName))
  }

  getPos.NAs = function(mat) {
    tmp = apply_mat_filter(mat)
    ind = which(is.na(tmp), arr.ind = TRUE)
    Pos = ind
    Pos[, 1] = ind[, 2]
    Pos[, 2] = -ind[, 1] + 1 + n
    return(Pos)
  }

  if (is.null(col)) {
    col = corrplot::COL2("RdBu", 200)
  }

  expand_expression = function(s) {
    ifelse(grepl("^[:=$]", s), parse(text = substring(s,
                                                      2)), s)
  }

  assign.color = function(dat = DAT, color = col) {
    newcorr = dat
    newcorr[newcorr <= 0] = 0
    newcorr[newcorr >= 1] = 1 - 0.0000000000000001
    color[floor(newcorr * length(color)) + 1]
  }

  pie.dat = function(theta, length = 100) {
    k = seq(pi/2, pi/2 - theta, length = theta * length *
              abs(theta)/pi)
    x = c(0, cos(k)/2, 0)
    y = c(0, sin(k)/2, 0)
    cbind(rbind(x, y), c(NA, NA))
  }

  isFALSE = function(x) identical(x, FALSE)
  isTRUE = function(x) identical(x, TRUE)

  ### COLOURS ####

  col.lim = c(-1, 1)
  intercept = 0
  zoom = 1
  col.lim2 = (intercept + col.lim) * zoom
  int = intercept * zoom

  if (is.null(col)) {
    col = COL2("RdBu", 200)
  }

  ### MATRIX SPECIFICATIONS ####

  n = nrow(corr)
  m = ncol(corr)
  min.nm = min(n, m)
  ord = 1:min.nm

  if (is.null(rownames(corr))) {
    rownames(corr) = 1:n
  }

  if (is.null(colnames(corr))) {
    colnames(corr) = 1:m
  }

  testTemp = getPos.Dat(corr)
  Pos = getPos.Dat(corr)[[1]]
  PosName = getPos.Dat(corr)[[3]]

  AllCoords = rbind(Pos)
  n2 = max(AllCoords[, 2])
  n1 = min(AllCoords[, 2])
  nn = n2 - n1
  m2 = max(AllCoords[, 1])
  m1 = min(AllCoords[, 1])
  mm = max(1, m2 - m1)

  newrownames = sapply(rownames(corr)[(n + 1 - n2):(n + 1 - n1)], expand_expression)
  newcolnames = sapply(colnames(corr)[m1:m2], expand_expression)

  DAT = getPos.Dat(corr)[[2]]
  len.DAT = length(DAT)
  rm(expand_expression)

  DATprop = getPos.Dat(corrProp)[[2]]

  ### PLOTTING SPECIFICATIONS ####

  col.fill = assign.color()

  if (isFALSE(outline)) {
    col.border = col.fill
  } else if (isTRUE(outline)) {
    col.border = "black"
  } else if (is.character(outline)) {
    col.border = outline
  }  else {
    stop("Unsupported value type for parameter outline")
  }

  oldpar = par(mar = mar, bg = par()$bg)
  on.exit(par(oldpar), add = TRUE)


  ### PLOT ####

  if (!add) {
    plot.new()
    xlabwidth = max(strwidth(newrownames, cex = tl.cex))
    ylabwidth = max(strwidth(newcolnames, cex = tl.cex))
    laboffset = strwidth("W", cex = tl.cex) * tl.offset
    for (i in 1:50) {
      xlim = c(m1 - 0.5 - laboffset - xlabwidth * (grepl("l",
                                                         tl.pos) | grepl("d", tl.pos)), m2 + 0.5 + mm *
                 cl.ratio * (cl.pos == "r") + xlabwidth * abs(cos(tl.srt *
                                                                    pi/180)) * grepl("d", tl.pos))
      ylim = c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b") -
                 laboffset, n2 + 0.5 + laboffset + ylabwidth *
                 abs(sin(tl.srt * pi/180)) * grepl("t", tl.pos) +
                 ylabwidth * abs(sin(tl.srt * pi/180)) * (type ==
                                                            "lower") * grepl("d", tl.pos))
      plot.window(xlim, ylim, asp = 1, xaxs = "i", yaxs = "i")
      x.tmp = max(strwidth(newrownames, cex = tl.cex))
      y.tmp = max(strwidth(newcolnames, cex = tl.cex))
      laboffset.tmp = strwidth("W", cex = tl.cex) * tl.offset
      if (max(x.tmp - xlabwidth, y.tmp - ylabwidth, laboffset.tmp -
              laboffset) < 0.001) {
        break
      }
      xlabwidth = x.tmp
      ylabwidth = y.tmp
      laboffset = laboffset.tmp
      if (i == 50) {
        warning(c("Not been able to calculate text margin, ",
                  "please try again with a clean new empty window using ",
                  "{plot.new(); dev.off()} or reduce tl.cex"))
      }
    }
    if (.Platform$OS.type == "windows") {
      grDevices::windows.options(width = 7, height = 7 *
                                   diff(ylim)/diff(xlim))
    }

    xlim = xlim + diff(xlim) * 0.01 * c(-1, 1)
    ylim = ylim + diff(ylim) * 0.01 * c(-1, 1)

    plot.window(xlim = xlim, ylim = ylim,
                xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  }

  laboffset = strwidth("W", cex = tl.cex) * tl.offset

  symbols(Pos, add = TRUE, inches = FALSE, rectangles = matrix(1, len.DAT, 2), bg = bg, fg = bg)

  # Circles and pies

  PIE.dat = lapply(DATprop * 2 * pi, pie.dat)
  len.pie = unlist(lapply(PIE.dat, length))/2
  PIE.dat2 = matrix(unlist(PIE.dat), ncol = 2, byrow = TRUE)
  PIE.dat2 <- PIE.dat2 * rep(DAT, len.pie)
  PIE.dat2 <- PIE.dat2 * 0.85
  PIE.dat2 = PIE.dat2 + Pos[rep(1:length(DAT), len.pie), ]

  # symbols(Pos, add = TRUE, inches = FALSE, circles = rep(0.5, len.DAT) * 0.85, fg = col.border)
  symbols(Pos, add = TRUE, inches = FALSE, circles = 0.85 * abs(DAT)/2, fg = col.fill, bg = bg, lwd = 3)

  polygon(PIE.dat2, border = col.fill, col = col.fill)

  draw_grid(AllCoords, "black")


  # Legend
  if (cl.pos != "n") {
    colRange = assign.color(dat = col.lim2)
    ind1 = which(col == colRange[1])
    ind2 = which(col == colRange[2])
    colbar = col[ind1:ind2]
    if (is.null(cl.length)) {
      cl.length = ifelse(length(colbar) > 20, 11, length(colbar) +
                           1)
    }
    labels = seq(col.lim[1], col.lim[2], length = cl.length)
    if (cl.pos == "r") {
      vertical = TRUE
      xlim = c(m2 + 0.5 + mm * 0.02, m2 + 0.5 + mm * cl.ratio)
      ylim = c(n1 - 0.5, n2 + 0.5)
    }
    if (cl.pos == "b") {
      vertical = FALSE
      xlim = c(m1 - 0.5, m2 + 0.5)
      ylim = c(n1 - 0.5 - nn * cl.ratio, n1 - 0.5 - nn *
                 0.02)
    }
    corrplot::colorlegend(colbar = colbar, labels = round(labels, 2),
                          offset = cl.offset, ratio.colbar = 0.3, cex = cl.cex,
                          xlim = xlim, ylim = ylim, vertical = vertical, align = cl.align.text)
  }

  # labels
  if (tl.pos != "n") {
    pos.xlabel = cbind(m1:m2, n2 + 0.5 + laboffset)
    pos.ylabel = cbind(m1 - 0.5, n2:n1)
    if (tl.pos == "td") {
      if (type != "upper") {
        stop("type should be 'upper' if tl.pos is 'dt'.")
      }
      pos.ylabel = cbind(m1:(m1 + nn) - 0.5, n2:n1)
    }
    if (tl.pos == "ld") {
      if (type != "lower") {
        stop("type should be 'lower' if tl.pos is 'ld'.")
      }
      pos.xlabel = cbind(m1:m2, n2:(n2 - mm) + 0.5 + laboffset)
    }
    if (tl.pos == "d") {
      pos.ylabel = cbind(m1:(m1 + nn) - 0.5, n2:n1)
      pos.ylabel = pos.ylabel[1:min(n, m), ]
      symbols(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], add = TRUE,
              bg = bg, fg = addgrid.col, inches = FALSE, squares = rep(1,
                                                                       length(pos.ylabel[, 1])))
      text(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], newcolnames[1:min(n,
                                                                     m)], col = tl.col, cex = tl.cex)
    }
    else {
      if (tl.pos != "l") {
        text(pos.xlabel[, 1], pos.xlabel[, 2], newcolnames,
             srt = tl.srt, adj = ifelse(tl.srt == 0, c(0.5,
                                                       0), c(0, 0)), col = tl.col, cex = tl.cex,
             offset = tl.offset)
      }
      text(pos.ylabel[, 1], pos.ylabel[, 2], newrownames,
           col = tl.col, cex = tl.cex, pos = 2, offset = tl.offset)
    }
  }

  title(title)
}
