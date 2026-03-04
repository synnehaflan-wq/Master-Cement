library(dplyr)

# Parametere

Base_Parameter <- list(
  ## Etterspørsel
  eps = 0.5,            # ε: priselastisitet i etterspørselen (ε > 0)
  A   = 1e6,            # A: skaleringsparameter for etterspørselen (nivå/markedsstørrelse)
  
  ## Karbonpriser
  P_CO2  = 75,          # P_CO2: karbonpris i EU ETS (€/tCO2)
  P_home = 0,           # P_home: karbonpris i eksportlandet/ROW (€/tCO2)
   
  ## Utslippsintensitet
  I_eu  = 0.60,         # I_eu: utslippsintensitet EU (tCO2 per tonn sement)
  I_no  = 0.60,         # I_no: utslippsintensitet Norge (tCO2 per tonn sement)
  I_row = 0.60,         # I_row: utslippsintensitet ROW (tCO2 per tonn sement)
  
  ## CCS (karbonfangst og -lagring)
  alpha_eu = 0.00,      # α_eu: fangstandel EU (0–1) = andel utslipp fanget med CCS
  alpha_no = 0.42,      # α_no: fangstandel Norge (0–1)
  C_CCS_eu = 80,        # C_CCS_eu: CCS-kostnad EU (€/tonn sement) – fangst+transport+lagring
  C_CCS_no = 90,        # C_CCS_no: CCS-kostnad Norge (€/tonn sement)
  
  ## Gratiskvoter (output-basert tildeling)
  beta = 0.20,  # β: gratistildeling (tCO2 i gratiskvoter per tonn sement) -> verdi: P_CO2*β (€/tonn)
  
 gamma_cbam = 0,   # 0 = ingen CBAM i basisåret, 1 = CBAM på # vet ikke om dette funker men variabel for å 1 ha CBAM på eller 0 ikke CBAM
  
  ## Kostnadsparametere i marginalkostnad: MC_r(x) = C0_r + C1_r*x + (policyledd)
  C0_eu  = 50,          # C0_eu: basekostnad EU (€/tonn) – konstantledd i MC
  C1_eu  = 1.1 ,   # C1_eu: helning EU (€/tonn^2) – økende MC når produksjon øker
  
  C0_no  = 55,          # C0_no: basekostnad Norge (€/tonn)
  C1_no  = 1.1  ,     # C1_no: helning Norge (€/tonn^2)
  
  C0_row = 45,          # C0_row: basekostnad ROW (€/tonn)
  C1_row = 1.1        # C1_row: helning ROW (€/tonn^2)
)


#____________________________________

# Likningene

## Marginal kostnadsfunksjoner MC_x
MC_eu <- function(x, p){ #EU sin MC
  p$C0_eu +
    p$C1_eu * x +
    p$P_CO2 * p$I_eu * (1 - p$alpha_eu) +
    p$alpha_eu * p$C_CCS_eu -
    p$P_CO2 * p$beta
}

MC_no <- function(x, p){ #Norge sin MC
  p$C0_no +
    p$C1_no * x +
    p$P_CO2 * p$I_no * (1 - p$alpha_no) +
    p$alpha_no * p$C_CCS_no -
    p$P_CO2 * p$beta
}

MC_row <- function(x, p){  #resten av verden sin MC
  p$C0_row +
    p$C1_row * x +
    p$gamma_cbam * p$P_CO2 * p$I_row  #fikse dette imorgen  dette kan fungere som referanse senario men den vil ikke ta med kvotepris ettersom P_home,men midlertidig løsning for å få frem 
}  

## Etterspørsels likning
Qd <- function(P, p){
  p$A * P^(-p$eps)
}


## Tilbud 
# Finner produksjon x slik at:
# P = MC(x)
# Fullkommen konkurranse -> produsenten produserer til pris = marginalkostnad

supply_from_MC <- function(P, MC_fun, p){
  
  if (MC_fun(0, p) >= P) return(0)
  # Hvis prisen er lavere enn marginalkostnad ved null produksjon,
  # produserer bedriften ingenting
  
  f <- function(x) MC_fun(x, p) - P
  # omskriver ligningen MC(x) = P til MC(x) - P = 0
  
  upper <- 1
  while (f(upper) < 0) upper <- upper * 2
  # øker øvre grense til MC overstiger P
  
  uniroot(f, lower = 0, upper = upper)$root
  # Numerisk løsning på MC(x) = P
}

## Markeds klarering
# Excess demand = Etterspørsel - Tilbud
# Likevekt når dette = 0

excess_demand <- function(P, p){
  
  x_eu  <- supply_from_MC(P, MC_eu, p)
  x_no  <- supply_from_MC(P, MC_no, p)
  x_row <- supply_from_MC(P, MC_row, p)
  
  Qd(P, p) - (x_eu + x_no + x_row)
}

## løse likevekt
# Finner markedspris P slik at:
# Q_D(P) = x_eu + x_no + x_row

solve_equilibrium <- function(p){
  
  P_star <- uniroot(function(P) excess_demand(P, p),
                    lower = 1e-6, upper = 5000)$root
  
  # Når P* er funnet, finner vi tilhørende produksjon
  
  x_eu  <- supply_from_MC(P_star, MC_eu, p)
  x_no  <- supply_from_MC(P_star, MC_no, p)
  x_row <- supply_from_MC(P_star, MC_row, p)
  
  list(P = P_star,
       x_eu = x_eu,
       x_no = x_no,
       x_row = x_row,
       Q = Qd(P_star, p))
}


# ============================================================
# KALIBRERING TIL BASISÅR 2024 (Reference-policy i 2024)
# Enheter: Q og x i Mt (millioner tonn), P i €/tonn
# Plasseres rett etter solve_equilibrium()
# ============================================================

calib_2024 <- list(
  P_target     = 160,        # €/tonn  (BYTT til 2024-pris)
  Q_target     = 160,        # Mt      (BYTT til 2024-forbruk/etterspørsel)
  x_row_target = 11.369,     # Mt      (BYTT til 2024-import ROW)
  x_no_target  = 1.132       # Mt      (BYTT til 2024-Norge forbruk)
)

# EU-leveranse = residual (må være >= 0)
calib_2024$x_eu_target <- calib_2024$Q_target -
  calib_2024$x_row_target - calib_2024$x_no_target
stopifnot(calib_2024$x_eu_target >= 0)

# (1) Kalibrer A slik at Q(P_target) = Q_target
# Q = A * P^(-eps)  =>  A = Q * P^(eps)
Base_Parameter$A <- calib_2024$Q_target * (calib_2024$P_target ^ Base_Parameter$eps)

# (2) Policy-ledd i basisåret (bruker Base_Parameter sine 2024-policyverdier)
policy_eu_2024 <- Base_Parameter$P_CO2 * Base_Parameter$I_eu * (1 - Base_Parameter$alpha_eu) +
  Base_Parameter$alpha_eu * Base_Parameter$C_CCS_eu -
  Base_Parameter$P_CO2 * Base_Parameter$beta

policy_no_2024 <- Base_Parameter$P_CO2 * Base_Parameter$I_no * (1 - Base_Parameter$alpha_no) +
  Base_Parameter$alpha_no * Base_Parameter$C_CCS_no -
  Base_Parameter$P_CO2 * Base_Parameter$beta

# 2024: uten CBAM (ROW betaler bare hjemlandspris, ofte 0)
policy_row_2024 <- Base_Parameter$P_home * Base_Parameter$I_row

# (3) Direkte C0-kalibrering: P_target = C0 + C1*x_target + policy
Base_Parameter$C0_eu  <- calib_2024$P_target - Base_Parameter$C1_eu  * calib_2024$x_eu_target  - policy_eu_2024
Base_Parameter$C0_no  <- calib_2024$P_target - Base_Parameter$C1_no  * calib_2024$x_no_target  - policy_no_2024
Base_Parameter$C0_row <- calib_2024$P_target - Base_Parameter$C1_row * calib_2024$x_row_target - policy_row_2024

# (4) Sjekk basisåret
eq_2024_check <- solve_equilibrium(Base_Parameter)

cat("\n===== SJEKK BASISÅR 2024 (kalibrert) =====\n")
cat("P (modell):", round(eq_2024_check$P, 4), " | P_target:", calib_2024$P_target, "\n")
cat("Q (modell):", round(eq_2024_check$Q, 4), " | Q_target:", calib_2024$Q_target, "\n")
cat("x_eu (modell):", round(eq_2024_check$x_eu, 4), " | x_eu_target:", calib_2024$x_eu_target, "\n")
cat("x_no (modell):", round(eq_2024_check$x_no, 4), " | x_no_target:", calib_2024$x_no_target, "\n")
cat("x_row (modell):", round(eq_2024_check$x_row, 4), " | x_row_target:", calib_2024$x_row_target, "\n")
cat("=========================================\n\n")
#__________________________________________________________

# Senarioer

## BAU (buisnis as usual) 
# Scenario uten noen form for regulering:
# - Ingen EU ETS (P_CO2 = 0)
# - Ingen gratiskvoter (beta = 0)
# - Ingen CBAM (P_home irrelevant når P_CO2 = 0)
# - Ingen CCS (alpha = 0)
#
# Økonomisk tolkning:
# Dette er et frimarked uten klimaregulering.
# Alle produsenter møter kun produksjonskostnader.
# Karbon har ingen pris.
# Prisen bestemmes kun av C0, C1 og etterspørsel.

Scenario_BaU <- Base_Parameter

Scenario_BaU$P_CO2 <- 0        # Ingen karbonpris i EU
Scenario_BaU$P_home <- 0       # Ingen karbonpris i ROW
Scenario_BaU$beta <- 0         # Ingen gratiskvoter
Scenario_BaU$alpha_eu <- 0     # Ingen CCS i EU
Scenario_BaU$alpha_no <- 0     # Ingen CCS i Norge



## Referanse senario
# - EU ETS aktiv (P_CO2 > 0)
# - Gratiskvoter eksisterer (beta > 0)
# - Ingen CBAM (ROW betaler ikke EU-karbonpris)
# - Norge har CCS
# - EU har ikke CCS
#
# Økonomisk tolkning:
# EU og Norge betaler karbonpris,
# men får delvis kompensasjon gjennom gratiskvoter.
# Norge har lavere effektiv utslippsintensitet grunnet CCS.
# ROW har konkurransefordel fordi de ikke betaler karbonpris.

Scenario_Reference <- Base_Parameter

Scenario_Reference$P_CO2 <- 75        # EU ETS aktiv
Scenario_Reference$beta <- 0.2        # Gratiskvoter (eksempelverdi)
Scenario_Reference$alpha_eu <- 0      # Ingen CCS i EU
Scenario_Reference$alpha_no <- 0.42   # CCS i Norge
Scenario_Reference$P_home <- 0        # Ingen karbonpris i ROW


#_____________________________________________________________________
# SIMULATE_PATH(): tidsserie 2025–2035 (serie av statiske likevekter)
# - Vi simulerer ett scenario år for år
# - Hvert år: vi oppdaterer policy-parametere (typisk beta, evt. alpha_eu)
# - Vi løser likevekt for hvert år og lagrer resultater
#
# Dette gir en tabell dere kan plotte (pris, produksjon, import, osv.)

simulate_path <- function(par_base,
                          scenario_name,
                          year_start = 2025,
                          year_end = 2035,
                          beta_start = NULL,
                          beta_end = NULL,
                          alpha_eu_path = NULL,
                          cbam_on_from_2026 = FALSE){
  
  years <- year_start:year_end
  n <- length(years)
  
  
  # 1) Lag beta-bane hvis ønskelig(mulig nødvendig)
  
  # Hvis beta_start/beta_end er oppgitt: lag lineær utfasing.
  # Hvis ikke: behold par_base$beta konstant i alle år.
  if (!is.null(beta_start) && !is.null(beta_end)) {
    beta_vec <- seq(beta_start, beta_end, length.out = n)
  } else {
    beta_vec <- rep(par_base$beta, n)
  }
  
  
  # 2) Lag alpha_eu-bane hvis ønskelig( tror kanskje nødvendig for å få till men usikker)
  
  # Hvis alpha_eu_path ikke er gitt: hold alpha_eu konstant.
  if (is.null(alpha_eu_path)) {
    alpha_eu_vec <- rep(par_base$alpha_eu, n)
  } else {
    if (length(alpha_eu_path) != n) stop("alpha_eu_path må ha samme lengde som antall år.")
    alpha_eu_vec <- alpha_eu_path
  }
  
  
  # 3) Simuler år-for-år
  
  out <- lapply(seq_along(years), function(i){
    
    # kopi av parameterne for dette året
    p <- par_base
    
    # oppdater policy-variabler
    p$beta <- beta_vec[i]
    p$alpha_eu <- alpha_eu_vec[i]
    
    # CBAM: enkel av/på (ikke gradvis vi må kanskje legge inn et gamma ledd)
    # - 2025: CBAM av (gjør (P_CO2 - P_home)=0)
    # - 2026+: CBAM på (ROW møter EU-pris via MC_row)
    if (cbam_on_from_2026) {
      if (years[i] <= 2025) {
        p$P_home <- p$P_CO2  # gir (P_CO2 - P_home)=0 -> ingen CBAM
      } else {
        p$P_home <- 0        # full CBAM fra 2026
      }
    }
    
    # løs likevekt
    eq <- solve_equilibrium(p)
    
    # lagre resultater
    data.frame(
      scenario = scenario_name,
      year = years[i],
      beta = p$beta,
      alpha_eu = p$alpha_eu,
      P = eq$P,
      Q = eq$Q,
      x_eu = eq$x_eu,
      x_no = eq$x_no,
      x_row = eq$x_row,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}


# BaU: alt konstant (ingen bane nødvendig)
path_BaU <- simulate_path(Scenario_BaU, "BaU", year_start = 2025, year_end = 2035)

# Referanse: beta fases ut (beta se litt mer på)
path_REF <- simulate_path(
  Scenario_Reference, "Reference",
  year_start = 2025, year_end = 2035,
  beta_start = 0.2, beta_end = 0,
  cbam_on_from_2026 = FALSE   # sett TRUE hvis vil ha CBAM "på" fra 2026 (ikke gradvis kan hende vi må legge inn ledd for å få i gradevis infasing)
)