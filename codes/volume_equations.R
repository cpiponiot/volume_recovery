#######################################################################################################################
##################################             VOLUME EQUATIONS                      ##################################
#######################################################################################################################

load("C:/Users/camille.piponiot/Google Drive/volume/gfbi/data/pars_max.Rdata")


volume = function(t,
                  ag = 4.4,
                  am = 1.4,
                  bg = 0.0116,
                  bm = 0.0109,
                  th = 0.015) {
  ag / th * (1 - (th * exp(-bg * t) - bg * exp(-th * t)) / (th - bg)) - am /
    th * (1 - (th * exp(-bm * t) - bm * exp(-th * t)) / (th - bm))
}

inv_volume = function(v,
                      ag = 4.4,
                      am = 1.4,
                      bg = 0.0116,
                      bm = 0.0109,
                      th = 0.015,
                      tmax = 500) {
  if (v < (ag - am) / th) {
    uniroot((function (x)
      volume(x, ag, am, bg, bm, th) - v), lower = 0, upper = tmax)[1]
  } else
    NA
}

dVG = function(t,
               ag = 4.4,
               am = 1.4,
               bg = 0.0116,
               bm = 0.0109,
               th = 0.015) {
  ag * bg / (th - bg) * (exp(-bg * t) - exp(-th * t)) + am * (1 - (th * exp(-bm *
                                                                              t) - bm * exp(-th * t)) / (th - bm))
}

dVM = function(t, am = pars_max$aM, bm = pars_max$bM) {
  am * (1 - exp(-bm * t))
}

cum_dVG = function(t,
                   t0 = 150,
                   ag = 4.4,
                   am = 1.4,
                   bg = 0.0116,
                   bm = 0.0109,
                   th = 0.015) {
  ag * bg / (th - bg) * ((exp(-bg * t0) - exp(-bg * t)) / bg - (exp(-th *
                                                                      t0) - exp(-th * t)) / th) +
    am * ((t - t0) - (th / bm * (exp(-bm * t0) - exp(-bm * t)) - bm / th *
                        (exp(-th * t0) - exp(-th * t))) / (th - bm))
}

cum_dVM = function(t,
                   t0 = 150,
                   am = 1.4,
                   bm = 0.0109) {
  am * ((t - t0) - (exp(-bm * t0) - exp(-bm * t)) / bm)
}
