vol_equations <- function(dbh, site, E){
  site = unique(site)
  vol = rep(0, length(dbh))
  if (site == "ita"){ vol <- 0.001602*dbh^1.9}
  else {if (site %in% c("ira", "tbc", "cum", "chb","bsl")){vol <- 0.000308*(dbh)^2.1988 }
    else {if (site %in% c("eco", "pet","prg") ) {vol <- 10^(-2.96 + 1.93 *log10(dbh))}
      else {if (site =="tpj"){vol <- exp(-7.62812 + 2.18090 * log(dbh))}
        else {if (site =="jar"){vol <- (-0.367921 + 0.0013446*(dbh)^2)}
          else {if (site=="prc"){vol <- -0.035829+0.00087634*(dbh)^2 }
            else { vol <- 0.65*exp(0.893 - E + 0.760*log(dbh) - 0.0340*(log(dbh))^2)*((dbh/200)^2*pi)}
          } 
        } 
      } 
    } 
  } 
  return(vol)
}