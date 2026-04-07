library(ggplot2)
library(tidyr)
library(dplyr)


# Parametere

Base_Parameter <- list(
  ## Etterspørsel
  eps = 0.3,            # ε: priselastisitet i etterspørselen (ε > 0)
  A   = 1e6,            # A: skaleringsparameter for etterspørselen (nivå/markedsstørrelse)
  
  ## Karbonpriser
  P_CO2  = 65,          # P_CO2: karbonpris i EU ETS (€/tCO2)
  P_home = 0,           # P_home: karbonpris i eksportlandet/ROW (€/tCO2)
   
  ## Utslippsintensitet
  I_eu  = 0.72,         # I_eu: utslippsintensitet EU (tCO2 per tonn sement)
  I_no  = 0.72,         # I_no: utslippsintensitet Norge (tCO2 per tonn sement)
  I_row = 0.72,         # I_row: utslippsintensitet ROW (tCO2 per tonn sement)
  
  ## CCS (karbonfangst og -lagring)
  alpha_eu = 0.00,      # α_eu: fangstandel EU (0–1) = andel utslipp fanget med CCS
  alpha_no = 0.0,      # α_no: fangstandel Norge (0–1)
  C_CCS_eu = 180,        # C_CCS_eu: CCS-kostnad EU (€/tonn sement) – fangst+transport+lagring
  C_CCS_no = 180,        # C_CCS_no: CCS-kostnad Norge (€/tonn sement)
  
  ## Gratiskvoter (output-basert tildeling)
  beta = 0.69,  # β: gratistildeling (tCO2 i gratiskvoter per tonn sement) -> verdi: P_CO2*β (€/tonn)
  
 gamma_cbam = 0,   # 0 = ingen CBAM i basisåret, 1 = CBAM på # vet ikke om dette funker men variabel for å 1 ha CBAM på eller 0 ikke CBAM
  
 ## subsidie på CCS 
 s_eu = 0.65 ,
 s_no = 0.65 ,
 
 ## parameter for hvor stor andel av eu som har ccs
 rho_eu = 0 ,

 
  ## Kostnadsparametere i marginalkostnad: MC_r(x) = C0_r + C1_r*x + (policyledd)
  C0_eu  = 35.55,          # C0_eu: basekostnad EU (€/tonn) – konstantledd i MC
  C1_eu  = 0.5 ,   # C1_eu: helning EU (€/tonn^2) – økende MC når produksjon øker
  
  C0_no  = 35.55,          # C0_no: basekostnad Norge (€/tonn)
  C1_no  = 2.5  ,     # C1_no: helning Norge (€/tonn^2)
  
  C0_row = 35.55,          # C0_row: basekostnad ROW (€/tonn)
  C1_row = 3.0        # C1_row: helning ROW (€/tonn^2)
)

#____________________________________

# Likningene

## Marginal kostnadsfunksjoner MC_x
MC_eu <- function(x, p){ # EU sin MC
  p$C0_eu +
    p$C1_eu * x +
    p$P_CO2 * p$I_eu * (1 - p$alpha_eu) +
    p$alpha_eu * p$rho_eu * (1 - p$s_eu) * p$C_CCS_eu -
    p$P_CO2 * p$beta
}

MC_no <- function(x, p){ # Norge sin MC
  p$C0_no +
    p$C1_no * x +
    p$P_CO2 * p$I_no * (1 - p$alpha_no) +
    p$alpha_no * (1 - p$s_no) * p$C_CCS_no -
    p$P_CO2 * p$beta
}

MC_row <- function(x, p){  # Resten av verden sin MC
  p$C0_row +
    p$C1_row * x +
    p$P_home * p$I_row +
    p$gamma_cbam * (p$P_CO2 - p$P_home) * p$I_row
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
  P_target     = 136.1,      # €/tonn
  Q_target     = 327.9,    # Mt
  x_row_target = 11.6,   # Mt
  x_no_target  = 1.8     # Mt
)

# EU-leveranse = residual
calib_2024$x_eu_target <- calib_2024$Q_target -
  calib_2024$x_row_target - calib_2024$x_no_target

stopifnot(calib_2024$x_eu_target >= 0)

# Sjekk at ingen produksjonsnivåer er null
stopifnot(calib_2024$x_eu_target  > 0)
stopifnot(calib_2024$x_no_target  > 0)
stopifnot(calib_2024$x_row_target > 0)

# ------------------------------------------------------------
# (1) Kalibrer A slik at Q(P_target) = Q_target
# Q = A * P^(-eps)  =>  A = Q * P^(eps)
# ------------------------------------------------------------

Base_Parameter$A <- calib_2024$Q_target * (calib_2024$P_target ^ Base_Parameter$eps)

# ------------------------------------------------------------
# (2) Policy-ledd i basisåret
# ------------------------------------------------------------

policy_eu_2024 <- Base_Parameter$P_CO2 * Base_Parameter$I_eu * (1 - Base_Parameter$alpha_eu) +
  Base_Parameter$alpha_eu * Base_Parameter$rho_eu * (1 - Base_Parameter$s_eu) * Base_Parameter$C_CCS_eu -
  Base_Parameter$P_CO2 * Base_Parameter$beta


policy_no_2024 <- Base_Parameter$P_CO2 * Base_Parameter$I_no * (1 - Base_Parameter$alpha_no) +
  Base_Parameter$alpha_no * (1 - Base_Parameter$s_no) * Base_Parameter$C_CCS_no -
  Base_Parameter$P_CO2 * Base_Parameter$beta

policy_row_2024 <- Base_Parameter$P_home * Base_Parameter$I_row +
  Base_Parameter$gamma_cbam * (Base_Parameter$P_CO2 - Base_Parameter$P_home) * Base_Parameter$I_row

# ------------------------------------------------------------
# (3) Direkte C1-kalibrering
# ------------------------------------------------------------
# Vi holder C0 fast og kalibrerer helningen C1 slik at
# marginalkostnaden treffer observert pris i basisåret:
#
# P_target = C0 + C1*x_target + policy
# => C1 = (P_target - C0 - policy) / x_target
# ------------------------------------------------------------

Base_Parameter$C1_eu  <- (calib_2024$P_target - Base_Parameter$C0_eu  - policy_eu_2024)  / calib_2024$x_eu_target
Base_Parameter$C1_no  <- (calib_2024$P_target - Base_Parameter$C0_no  - policy_no_2024)  / calib_2024$x_no_target
Base_Parameter$C1_row <- (calib_2024$P_target - Base_Parameter$C0_row - policy_row_2024) / calib_2024$x_row_target

# ------------------------------------------------------------
# (4) Sjekk basisåret
# ------------------------------------------------------------

eq_2024_check <- solve_equilibrium(Base_Parameter)

cat("\n===== SJEKK BASISÅR 2024 (kalibrert) =====\n")
cat("P (modell):", round(eq_2024_check$P, 4), " | P_target:", calib_2024$P_target, "\n")
cat("Q (modell):", round(eq_2024_check$Q, 4), " | Q_target:", calib_2024$Q_target, "\n")
cat("x_eu (modell):", round(eq_2024_check$x_eu, 4), " | x_eu_target:", calib_2024$x_eu_target, "\n")
cat("x_no (modell):", round(eq_2024_check$x_no, 4), " | x_no_target:", calib_2024$x_no_target, "\n")
cat("x_row (modell):", round(eq_2024_check$x_row, 4), " | x_row_target:", calib_2024$x_row_target, "\n")
cat("C1_eu:", round(Base_Parameter$C1_eu, 6), "\n")
cat("C1_no:", round(Base_Parameter$C1_no, 6), "\n")
cat("C1_row:", round(Base_Parameter$C1_row, 6), "\n")
cat("=========================================\n\n")



#==============================================
# Senarioer liste
#===============================================


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
Scenario_BaU$gamma_cbam <- 0


# ============================================================
# REFERANSESCENARIO 1
# ------------------------------------------------------------
# Markedssituasjon tilsvarende 2024
#
# - EU ETS aktiv
# - Gratiskvoter eksisterer
# - Ingen CBAM
# - Ingen CCS i EU
# - Ingen CCS i Norge
#
# Økonomisk tolkning:
# Produsenter i EU og Norge møter karbonpris,
# men mottar gratiskvoter som reduserer kostnaden.
# Produsenter i resten av verden betaler ingen karbonpris.
# ============================================================

Scenario_REF1 <- Base_Parameter

Scenario_REF1$P_CO2 <- 65
Scenario_REF1$beta <- 0.69
Scenario_REF1$alpha_eu <- 0
Scenario_REF1$alpha_no <- 0
Scenario_REF1$P_home <- 0
Scenario_REF1$gamma_cbam <- 0


# ============================================================
# REFERANSESCENARIO 2
# ------------------------------------------------------------
# Markedssituasjon med CCS i Norge
#
# - EU ETS aktiv
# - Gratiskvoter eksisterer
# - Ingen CBAM
# - Ingen CCS i EU
# - CCS i Norge
#
# Økonomisk tolkning:
# Norske produsenter reduserer utslipp gjennom CCS
# og får dermed lavere effektiv karbonkostnad.
# Alle andre rammebetingelser er de samme som i referansescenario 1.
# ============================================================

Scenario_REF2 <- Base_Parameter

Scenario_REF2$P_CO2 <- 65
Scenario_REF2$beta <- 0.69
Scenario_REF2$alpha_eu <- 0
Scenario_REF2$alpha_no <- 0.42
Scenario_REF2$P_home <- 0
Scenario_REF2$gamma_cbam <- 0

# ============================================================
# LØS LIKEVEKT FOR BASELINE-SCENARIOER
# ============================================================

eq_BaU  <- solve_equilibrium(Scenario_BaU)
eq_REF1 <- solve_equilibrium(Scenario_REF1)
eq_REF2 <- solve_equilibrium(Scenario_REF2)
#==============================================================


## tids senarioer(1-)
#==================================================

# Scenario 1:
# - EU ETS videreføres
# - Gratiskvoter fases helt ut
# - CBAM innføres på import til EU
# - Norge har CCS
# - EU har ikke CCS
#
# Økonomisk tolkning:
# EU-produsenter møter full karbonkostnad fordi gratiskvotene forsvinner.
# Norske produsenter møter også karbonpris, men får lavere effektive utslipp
# fordi CCS reduserer andelen utslipp som må betales for.
# Import fra ROW blir ilagt CBAM, slik at importørene møter en karbonkostnad
# tilsvarende EU ETS ved eksport til EU-markedet.

Scenario_1 <- Base_Parameter

# EU ETS videreføres
Scenario_1$P_CO2 <- 65
Scenario_1$beta <- 0# Ingen gratiskvoter i selve scenarioet (de fases ut i simuleringsbanen)
Scenario_1$alpha_eu <- 0# CCS: # EU har ikke CCS
Scenario_1$alpha_no <- 0.42# Norge har CCS
Scenario_1$P_home <- 0# Ingen hjemlig karbonpris i ROW
Scenario_1$gamma_cbam <- 0# CBAM-status i basis for scenarioet
# Vi setter 0 her og lar simulate_path styre når CBAM slås på



# ============================================================
# SCENARIO 2
# ------------------------------------------------------------
# Scenario 2:
# - EU ETS videreføres
# - Gratiskvoter fases helt ut
# - CBAM innføres på import til EU
# - Norge har CCS
# - EU innfører CCS i deler av sementindustrien
#
# Økonomisk tolkning:
# Alle produsenter i EU og Norge møter karbonprisen gjennom EU ETS.
# Samtidig fases gratiskvotene ut, slik at karbonkostnaden slår fullt inn.
# Import fra ROW blir ilagt CBAM, slik at utenlandske produsenter møter en
# karbonkostnad tilsvarende EU ETS ved eksport til EU.
#
# Forskjellen fra Scenario 1 er at deler av EU-industrien også tar i bruk CCS.
# Dette reduserer gjennomsnittlig utslippsintensitet i EU og dermed den
# effektive karbonkostnaden for EU-produsenter.
# ============================================================

# ============================================================
# SCENARIO 2A–2F
# ------------------------------------------------------------
# Scenario 2 analyseres som flere alternative konfigurasjoner
# av CCS i EU. Variantene skiller seg i:
# - rho_eu: andel EU-produsenter som investerer i CCS
# - s_eu: andel av CCS-kostnaden som dekkes av subsidie
#
# Felles antakelser:
# - EU ETS videreføres
# - Gratiskvoter fases ut
# - CBAM innføres fra 2026
# - Norge har CCS
# - Fangstgrad i EU settes separat i path-blokkene
# ============================================================

Scenario_2A <- Base_Parameter
Scenario_2A$P_CO2 <- 75
Scenario_2A$beta <- 0.69
Scenario_2A$alpha_eu <- 0
Scenario_2A$alpha_no <- 0.42
Scenario_2A$rho_eu <- 0.26
Scenario_2A$s_eu <- 0.11
Scenario_2A$P_home <- 0
Scenario_2A$gamma_cbam <- 0

Scenario_2B <- Base_Parameter
Scenario_2B$P_CO2 <- 75
Scenario_2B$beta <- 0.69
Scenario_2B$alpha_eu <- 0
Scenario_2B$alpha_no <- 0.42
Scenario_2B$rho_eu <- 0.44
Scenario_2B$s_eu <- 0.05
Scenario_2B$P_home <- 0
Scenario_2B$gamma_cbam <- 0

Scenario_2C <- Base_Parameter
Scenario_2C$P_CO2 <- 75
Scenario_2C$beta <- 0.69
Scenario_2C$alpha_eu <- 0
Scenario_2C$alpha_no <- 0.42
Scenario_2C$rho_eu <- 0.44
Scenario_2C$s_eu <- 0.4
Scenario_2C$P_home <- 0
Scenario_2C$gamma_cbam <- 0

Scenario_2D <- Base_Parameter
Scenario_2D$P_CO2 <- 75
Scenario_2D$beta <- 0.69
Scenario_2D$alpha_eu <- 0
Scenario_2D$alpha_no <- 0.42
Scenario_2D$rho_eu <- 0.78
Scenario_2D$s_eu <- 0.65
Scenario_2D$P_home <- 0
Scenario_2D$gamma_cbam <- 0

Scenario_2E <- Base_Parameter
Scenario_2E$P_CO2 <- 75
Scenario_2E$beta <- 0.69
Scenario_2E$alpha_eu <- 0
Scenario_2E$alpha_no <- 0.42
Scenario_2E$rho_eu <- 0.26
Scenario_2E$s_eu <- 0.05
Scenario_2E$P_home <- 0
Scenario_2E$gamma_cbam <- 0

Scenario_2F <- Base_Parameter
Scenario_2F$P_CO2 <- 75
Scenario_2F$beta <- 0.69
Scenario_2F$alpha_eu <- 0
Scenario_2F$alpha_no <- 0.42
Scenario_2F$rho_eu <- 0.8
Scenario_2F$s_eu <- 0.80
Scenario_2F$P_home <- 0
Scenario_2F$gamma_cbam <- 0


#_____________________________________________________________________
# SIMULATE_PATH(): tidsserie 2025–2035 (serie av statiske likevekter)
# --------------------------------------------------------------
# Vi simulerer ett scenario år for år.
#
# Hvert år gjør modellen følgende:
# 1. Oppdaterer policy-variabler (beta, alpha_eu, gamma_cbam)
# 2. Løser markedslikevekten gitt disse parameterne
# 3. Lagrer resultatene (pris, produksjon, import osv.)
#
# Resultatet blir en tabell med én rad per år som kan brukes til
# grafer eller videre analyse. viktigste og det vi fokuserer på blir forskjeld 2025 og 2035 

simulate_path <- function(par_base,
                          scenario_name,
                          year_start = 2025,
                          year_end = 2035,
                          beta_start = NULL,
                          beta_end = NULL,
                          alpha_eu_path = NULL,
                          alpha_no_path = NULL,
                          gamma_cbam_path = NULL){
  
  # Lager vektor med alle år i simuleringen
  years <- year_start:year_end
  
  # Antall år
  n <- length(years)
  
  
  # --------------------------------------------------------------
  # 1) Lag bane for gratiskvoter (beta)
  # --------------------------------------------------------------
  
  # Hvis start og sluttverdi er oppgitt lager vi lineær utfasing
  # Hvis ikke holder beta seg konstant i hele perioden
  
  if (!is.null(beta_start) && !is.null(beta_end)) {
    
    # lineær reduksjon av gratiskvoter
    beta_vec <- seq(beta_start, beta_end, length.out = n)
    
  } else {
    
    # beta konstant
    beta_vec <- rep(par_base$beta, n)
    
  }
  
  
  # --------------------------------------------------------------
  # 2) Lag bane for CCS i EU (alpha_eu)
  # --------------------------------------------------------------
  
  # Hvis ingen bane er spesifisert holder vi alpha_eu konstant
  
  if (is.null(alpha_eu_path)) {
    
    alpha_eu_vec <- rep(par_base$alpha_eu, n)
    
  } else {
    
    # sikkerhetssjekk: riktig lengde
    if (length(alpha_eu_path) != n) {
      stop("alpha_eu_path må ha samme lengde som antall år.")
    }
    
    alpha_eu_vec <- alpha_eu_path
    
  }
  
  # --------------------------------------------------------------
  # 2b) Lag bane for CCS i Norge (alpha_no)
  # --------------------------------------------------------------
  
  if (is.null(alpha_no_path)) {
    
    alpha_no_vec <- rep(par_base$alpha_no, n)
    
  } else {
    
    if (length(alpha_no_path) != n) {
      stop("alpha_no_path må ha samme lengde som antall år.")
    }
    
    alpha_no_vec <- alpha_no_path
    
  }
  # --------------------------------------------------------------
  # 3) Lag bane for CBAM (gamma_cbam)
  # --------------------------------------------------------------
  
  # gamma_cbam = 0  → ingen CBAM
  # gamma_cbam = 1  → CBAM aktiv
  
  if (is.null(gamma_cbam_path)) {
    
    # hvis ikke spesifisert: bruk samme verdi hele perioden
    gamma_cbam_vec <- rep(par_base$gamma_cbam, n)
    
  } else {
    
    # sikkerhetssjekk
    if (length(gamma_cbam_path) != n) {
      stop("gamma_cbam_path må ha samme lengde som antall år.")
    }
    
    gamma_cbam_vec <- gamma_cbam_path
    
  }
  
  
  # --------------------------------------------------------------
  # 4) Simuler markedet år for år
  # --------------------------------------------------------------
  
  out <- lapply(seq_along(years), function(i){
    
    # kopi av parameterlisten slik at vi kan endre verdier for dette året
    p <- par_base
    
    
    # --------------------------------------------------------------
    # Oppdater policy-parametere for året
    # --------------------------------------------------------------
    
    # gratiskvoter
    p$beta <- beta_vec[i]
    
    # CCS i EU
    p$alpha_eu <- alpha_eu_vec[i]
    
    #CCS i norge
    p$alpha_no <- alpha_no_vec[i]
    
    # CBAM-status
    p$gamma_cbam <- gamma_cbam_vec[i]
    
    
    # --------------------------------------------------------------
    # Løs markedslikevekt
    # --------------------------------------------------------------
    
    eq <- solve_equilibrium(p)
    
    
    # --------------------------------------------------------------
    # Lagre resultatene for dette året
    # --------------------------------------------------------------
    
    data.frame(
      
      scenario = scenario_name,
      year = years[i],
      
      beta = p$beta,
      alpha_eu = p$alpha_eu,
      gamma_cbam = p$gamma_cbam,
      
      P = eq$P,     # markedspris
      Q = eq$Q,     # total etterspørsel
      
      x_eu = eq$x_eu,    # EU-produksjon
      x_no = eq$x_no,    # norsk eksport
      x_row = eq$x_row,  # import fra resten av verden
      
      stringsAsFactors = FALSE
      
    )
    
  })
  
  
  # --------------------------------------------------------------
  # Slår sammen resultatene til én tabell
  # --------------------------------------------------------------
  
  do.call(rbind, out)
  
}
#=======================================================
# Senario tabellene
#-------------------------------------------------------
# BaU: alt konstant (ingen bane nødvendig)
path_BaU <- simulate_path(Scenario_BaU, "BaU", year_start = 2025, year_end = 2035)

# Referanse: beta fases ut (beta se litt mer på)
path_REF2 <- simulate_path(
  Scenario_REF2,
  "Reference 2",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end = 0.69,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),# CCS i Norge fra 2025
  gamma_cbam_path = rep(0, length(2025:2035))
)


path_REF1 <- simulate_path(
  Scenario_REF1,
  "Reference 1",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end   = 0.69,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0., length(2025:2035)), # CCS i Norge fra 2025
  gamma_cbam_path = rep(0, length(2025:2035))
)

# KJØRING AV SCENARIO 1
# ------------------------------------------------------------
# Antakelser:
# - Gratiskvoter er fjernet helt fra start, usikker på dette
# - CBAM innføres fra 2026
# - EU har ikke CCS
# - Norge har CCS
# ============================================================

path_S1 <- simulate_path(
  Scenario_1,
  "CBAM senario",
  year_start = 2025,
  year_end = 2035,
  
  # ingen gratiskvoter
  beta_start = 0,
  beta_end = 0,
  
  # EU har ikke CCS
  alpha_eu_path = rep(0, length(2025:2035)),
  
  # CBAM av i 2025, på fra 2026
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)


# ============================================================
# KJØRING AV SCENARIO 2A–2F
# ------------------------------------------------------------
# Her antar vi at fabrikkene i EU som faktisk har CCS,
# oppnår 50 % fangstgrad.
# alpha_eu tolkes derfor som fangstgrad på CCS-anleggene,
# mens rho_eu er andelen av EU-bedriftene som har CCS.
# ============================================================

path_S2A <- simulate_path(
  Scenario_2A,
  "2A",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end   = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S2B <- simulate_path(
  Scenario_2B,
  "2B",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end   = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S2C <- simulate_path(
  Scenario_2C,
  "2C",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end   = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S2D <- simulate_path(
  Scenario_2D,
  "2D",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end   = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S2E <- simulate_path(
  Scenario_2E,
  "2E",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end   = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S2F <- simulate_path(
  Scenario_2F,
  "2F",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end   = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)




# ============================================================
# SAMLE ALLE SCENARIOER I ÉN TABELL for analyse 
# ------------------------------------------------------------
# Denne tabellen samler alle simulerte tidsbaner slik at vi
# enkelt kan hente ut bestemte år til resultatdelen.
# ============================================================

all_paths <- bind_rows(
  path_BaU,
  path_REF1,
  path_REF2,
  path_S1,
  path_S2A,
  path_S2B,
  path_S2C,
  path_S2D,
  path_S2E,
  path_S2F
)


# TABELL FOR 2025 OG 2035
# ------------------------------------------------------------
# Henter ut startår og sluttår for hvert scenario.
# Dette gir en ryddig sammenligningstabell til drøftingen.
# ----------------------------------------------------------

results_2025_2035 <- all_paths %>%
  filter(year %in% c(2025, 2035)) %>%
  arrange(scenario, year) %>%
  mutate(
    P = round(P, 2),
    Q = round(Q, 2),
    x_eu = round(x_eu, 2),
    x_no = round(x_no, 2),
    x_row = round(x_row, 2)
  )

print(results_2025_2035)


# -------------------------------------------------------------
# ENDRING FRA 2025 TIL 2035
# ------------------------------------------------------------
# Denne tabellen viser hvor mye pris, etterspørsel og tilbud
# endrer seg mellom startåret og sluttåret i hvert scenario.
# ------------------------------------------------------------

changes_2025_2035 <- all_paths %>%
  filter(year %in% c(2025, 2035)) %>%
  select(scenario, year, P, Q, x_eu, x_no, x_row) %>%
  tidyr::pivot_wider(
    names_from = year,
    values_from = c(P, Q, x_eu, x_no, x_row)
  ) %>%
  mutate(
    dP     = round(P_2035 - P_2025, 2),
    dQ     = round(Q_2035 - Q_2025, 2),
    dx_eu  = round(x_eu_2035 - x_eu_2025, 2),
    dx_no  = round(x_no_2035 - x_no_2025, 2),
    dx_row = round(x_row_2035 - x_row_2025, 2)
  )

print(changes_2025_2035)


#tabell på ulike nivå av rho og s_eu 
#=====================================
scenario2_config <- data.frame(
  scenario = c("Scenario 2A", "Scenario 2B", "Scenario 2C",
               "Scenario 2D", "Scenario 2E", "Scenario 2F"),
  rho_eu = c(0.26, 0.44, 0.44, 0.78, 0.26, 0.80),
  s_eu   = c(0.11, 0.05, 0.4, 0.65, 0.05, 0.80)
)

print(scenario2_config)




# ============================================================
# GRAFDATA – SCENARIO 2A–2D SAMMENLIGNET MED REFERENCE 1
# ------------------------------------------------------------
# Vi bruker 2035-tallene og beregner prosentvis endring i
# pris (P) og etterspørsel (Q) relativt til Reference 1.
# ============================================================

plotdata_S2_vs_REF1 <- bind_rows(
  path_REF1,
  path_S2A,
  path_S2B,
  path_S2C,
  path_S2D
) %>%
  filter(year == 2035) %>%
  select(scenario, P, Q) %>%
  pivot_longer(
    cols = c(P, Q),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(variable) %>%
  mutate(
    ref_value = value[scenario == "Reference 1"],
    pct_change = 100 * (value - ref_value) / ref_value
  ) %>%
  ungroup() %>%
  filter(scenario != "Reference 1") %>%
  mutate(
    variable = recode(variable,
                      P = "Pris",
                      Q = "Etterspørsel"
    )
  )

print(plotdata_S2_vs_REF1)

#graf___________________________-

ggplot(plotdata_S2_vs_REF1, aes(x = scenario, y = pct_change, fill = variable)) +
  geom_col(
    position = position_dodge(width = 0.6),
    width = 0.5
  ) +
  geom_text(
    aes(
      label = paste0(round(pct_change, 1), "%"),
      y = ifelse(pct_change >= 0, pct_change + 1, pct_change - 1)
    ),
    position = position_dodge(width = 0.6),
    size = 4
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(
    data = distinct(plotdata_S2_vs_REF1, scenario),
    aes(x = scenario, y = 0, label = scenario),
    inherit.aes = FALSE,
    vjust = -0.8,
    size = 4
  ) +
  labs(
    title = "Scenario 2A–2D: endring i pris og etterspørsel i 2035 relativt til Reference 1",
    x = NULL,
    y = "Prosentvis endring fra Reference 1 (%)",
    fill = NULL
  ) +
  scale_y_continuous(
    breaks = seq(-15, 45, by = 5)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )
# ============================================================
# GRAFDATA – ABSOLUTT ENDRING MOT REFERENCE 1
# ============================================================

plotdata_S2_vs_REF1_abs <- bind_rows(
  path_REF1,
  path_S2A,
  path_S2B,
  path_S2C,
  path_S2D
) %>%
  filter(year == 2035) %>%
  select(scenario, P, Q) %>%
  pivot_longer(
    cols = c(P, Q),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(variable) %>%
  mutate(
    ref_value = value[scenario == "Reference 1"],
    abs_change = value - ref_value
  ) %>%
  ungroup() %>%
  filter(scenario != "Reference 1") %>%
  mutate(
    variable = recode(variable,
                      P = "Pris",
                      Q = "Etterspørsel"
    )
  )

ggplot(plotdata_S2_vs_REF1_abs, aes(x = scenario, y = abs_change, fill = variable)) +
  geom_col(position = "dodge") +
  labs(
    title = "Scenario 2A–2D: absolutt endring i pris og etterspørsel i 2035 relativt til Reference 1",
    x = NULL,
    y = "Absolutt endring fra Reference 1",
    fill = NULL
  ) +
  theme_minimal()
# ============================================================
# RUN_SENSITIVITY()
# ------------------------------------------------------------
# This function runs one sensitivity case for a chosen scenario.
# It updates selected parameter values, simulates the scenario,
# and returns the results.
# ============================================================

run_sensitivity <- function(par_scenario,
                            scenario_name,
                            case_name,
                            year_start = 2025,
                            year_end = 2035,
                            beta_start = NULL,
                            beta_end = NULL,
                            alpha_eu_path = NULL,
                            alpha_no_path = NULL,
                            gamma_cbam_path = NULL,
                            P_CO2_new = NULL,
                            eps_new = NULL,
                            C_CCS_eu_new = NULL,
                            C_CCS_no_new = NULL) {
  
  # Copy scenario
  p <- par_scenario
  
  # Update selected parameters if provided
  if (!is.null(P_CO2_new))    p$P_CO2 <- P_CO2_new
  if (!is.null(eps_new))      p$eps <- eps_new
  if (!is.null(C_CCS_eu_new)) p$C_CCS_eu <- C_CCS_eu_new
  if (!is.null(C_CCS_no_new)) p$C_CCS_no <- C_CCS_no_new
  
  # Keep subsidy tied to CCS cost if that is your assumption
  p$s_eu <- 0.65
  p$s_no <- 0.65
  
  # Run path
  out <- simulate_path(
    par_base = p,
    scenario_name = paste(scenario_name, "-", case_name),
    year_start = year_start,
    year_end = year_end,
    beta_start = beta_start,
    beta_end = beta_end,
    alpha_eu_path = alpha_eu_path,
    alpha_no_path = alpha_no_path,
    gamma_cbam_path = gamma_cbam_path
  )
  
  out$case <- case_name
  out
}








# ============================================================
# SENSITIVITY ANALYSIS – SCENARIO carbon cost (ETS pris)
# ------------------------------------------------------------
# We test how the results change when the carbon
# price is varied.

#- Scenario 1
# ============================================================

sens_S1_carbon_low <- run_sensitivity(
  par_scenario = Scenario_1,
  scenario_name = "Scenario 1",
  case_name = "Low carbon price",
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035))),
  P_CO2_new = 35
)

sens_S1_carbon_base <- run_sensitivity(
  par_scenario = Scenario_1,
  scenario_name = "Scenario 1",
  case_name = "Base carbon price",
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035))),
  P_CO2_new = 65
)

sens_S1_carbon_high <- run_sensitivity(
  par_scenario = Scenario_1,
  scenario_name = "Scenario 1",
  case_name = "High carbon price",
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035))),
  P_CO2_new = 145
)
sens_S1_carbon_highest <- run_sensitivity(
  par_scenario = Scenario_1,
  scenario_name = "Scenario 1",
  case_name = " Highest carbon price",
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035))),
  P_CO2_new = 200
)
sens_S1_carbon <- bind_rows(
  sens_S1_carbon_low,
  sens_S1_carbon_base,
  sens_S1_carbon_high,
  sens_S1_carbon_highest
)

sens_S1_carbon_2025_2035 <- sens_S1_carbon %>%
  filter(year %in% c(2025, 2035)) %>%
  select(case, year, P, Q, x_eu, x_no, x_row) %>%
  arrange(case, year) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  )

print(sens_S1_carbon_2025_2035)


# ============================================================
# SENSITIVITETSANALYSE – SCENARIO 2
# KARBONPRIS
# ------------------------------------------------------------
# Vi tester hvordan resultatene i Scenario 2 endres når
# karbonprisen varierer mellom lav, base og høy.
# ============================================================

sens_S2_carbon_low <- run_sensitivity(
  par_scenario = Scenario_2A,
  scenario_name = "Scenario 2",
  case_name = "Low carbon price",
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035))),
  P_CO2_new = 35
)

sens_S2_carbon_base <- run_sensitivity(
  par_scenario = Scenario_2A,
  scenario_name = "Scenario 2",
  case_name = "Base carbon price",
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035))),
  P_CO2_new = 75
)

sens_S2_carbon_high <- run_sensitivity(
  par_scenario = Scenario_2A,
  scenario_name = "Scenario 2",
  case_name = "High carbon price",
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035))),
  P_CO2_new = 145
)

sens_S2_carbon <- bind_rows(
  sens_S2_carbon_low,
  sens_S2_carbon_base,
  sens_S2_carbon_high
)

sens_S2_carbon_2025_2035 <- sens_S2_carbon %>%
  filter(year %in% c(2025, 2035)) %>%
  select(case, year, P, Q, x_eu, x_no, x_row) %>%
  arrange(case, year) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

print(sens_S2_carbon_2025_2035)


# ============================================================
# SENSITIVITETSANALYSE – CCS kostnad
# ------------------------------------------------------------
# Vi tester hvordan resultatene i Scenario 1 endres når
# CCS-kostnaden varierer mellom lav, base og høy.
#
# Antakelser:
# - Scenario 1 beholdes likt ellers
# - Kun CCS-kostnaden endres
# - Subsidien holdes lik 65 % av CCS-kostnaden
# ============================================================
# Senario 1
# --------
# Lav CCS-kostnad
# --------
Scenario_1_CCS_low <- Scenario_1
Scenario_1_CCS_low$C_CCS_eu <- 90
Scenario_1_CCS_low$C_CCS_no <- 90
Scenario_1_CCS_low$S_eu <- 0.65 * Scenario_1_CCS_low$C_CCS_eu
Scenario_1_CCS_low$S_no <- 0.65 * Scenario_1_CCS_low$C_CCS_no

path_S1_lowCCS <- simulate_path(
  Scenario_1_CCS_low,
  "Scenario 1",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S1_lowCCS$case <- "Low CCS cost"


# --------
# Basis CCS-kostnad
# --------
Scenario_1_CCS_base <- Scenario_1
Scenario_1_CCS_base$C_CCS_eu <- 180
Scenario_1_CCS_base$C_CCS_no <- 180
Scenario_1_CCS_base$S_eu <- 0.65 * Scenario_1_CCS_base$C_CCS_eu
Scenario_1_CCS_base$S_no <- 0.65 * Scenario_1_CCS_base$C_CCS_no

path_S1_baseCCS <- simulate_path(
  Scenario_1_CCS_base,
  "Scenario 1",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S1_baseCCS$case <- "Base CCS cost"


# --------
# Høy CCS-kostnad
# --------
Scenario_1_CCS_high <- Scenario_1
Scenario_1_CCS_high$C_CCS_eu <- 240
Scenario_1_CCS_high$C_CCS_no <- 240
Scenario_1_CCS_high$S_eu <- 0.65 * Scenario_1_CCS_high$C_CCS_eu
Scenario_1_CCS_high$S_no <- 0.65 * Scenario_1_CCS_high$C_CCS_no

path_S1_highCCS <- simulate_path(
  Scenario_1_CCS_high,
  "Scenario 1",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S1_highCCS$case <- "High CCS cost"


sens_S1_CCS <- bind_rows(
  path_S1_lowCCS,
  path_S1_baseCCS,
  path_S1_highCCS
)


# TABELL FOR 2025 OG 2035

sens_S1_CCS_2025_2035 <- sens_S1_CCS %>%
  filter(year %in% c(2025, 2035)) %>%
  select(case, year, P, Q, x_eu, x_no, x_row) %>%
  arrange(case, year) %>%
  mutate(
    P = round(P, 2),
    Q = round(Q, 2),
    x_eu = round(x_eu, 2),
    x_no = round(x_no, 2),
    x_row = round(x_row, 2)
  )

print(sens_S1_CCS_2025_2035)


# ============================================================
# SENSITIVITETSANALYSE – SCENARIO 2
# CCS-KOSTNAD
# ------------------------------------------------------------
# Vi tester hvordan resultatene i Scenario 2 endres når
# CCS-kostnaden varierer mellom lav, base og høy.
# ============================================================

Scenario_2_CCS_low <- Scenario_2A
Scenario_2_CCS_low$C_CCS_eu <- 90
Scenario_2_CCS_low$C_CCS_no <- 90

path_S2_lowCCS <- simulate_path(
  Scenario_2_CCS_low,
  "Scenario 2",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)
path_S2_lowCCS$case <- "Low CCS cost"


Scenario_2_CCS_base <- Scenario_2A
Scenario_2_CCS_base$C_CCS_eu <- 180
Scenario_2_CCS_base$C_CCS_no <- 180

path_S2_baseCCS <- simulate_path(
  Scenario_2_CCS_base,
  "Scenario 2",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)
path_S2_baseCCS$case <- "Base CCS cost"


Scenario_2_CCS_high <- Scenario_2A
Scenario_2_CCS_high$C_CCS_eu <- 240
Scenario_2_CCS_high$C_CCS_no <- 240

path_S2_highCCS <- simulate_path(
  Scenario_2_CCS_high,
  "Scenario 2",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)
path_S2_highCCS$case <- "High CCS cost"


sens_S2_CCS <- bind_rows(
  path_S2_lowCCS,
  path_S2_baseCCS,
  path_S2_highCCS
)

sens_S2_CCS_2025_2035 <- sens_S2_CCS %>%
  filter(year %in% c(2025, 2035)) %>%
  select(case, year, P, Q, x_eu, x_no, x_row) %>%
  arrange(case, year) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

print(sens_S2_CCS_2025_2035)

# ============================================================
# SENSITIVITETSANALYSE 
# ETTERSPØRSELSELASTISITET
# ------------------------------------------------------------
# Vi tester hvordan resultatene i Scenario 1 endres når
# etterspørselselastisiteten varierer mellom lav, base og høy.
#
# Viktig:
# Når eps endres, rekalibreres A slik at modellen fortsatt
# treffer observert etterspørsel i basisåret.
# ============================================================
#  – SCENARIO 1
# --------
# Lav elastisitet
# --------
Scenario_1_EPS_low <- Scenario_1
Scenario_1_EPS_low$eps <- 0.1
Scenario_1_EPS_low$A <- calib_2024$Q_target * (calib_2024$P_target ^ Scenario_1_EPS_low$eps)

path_S1_lowEPS <- simulate_path(
  Scenario_1_EPS_low,
  "Scenario 1",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S1_lowEPS$case <- "Low elasticity"


# --------
# Basis elastisitet
# --------
Scenario_1_EPS_base <- Scenario_1
Scenario_1_EPS_base$eps <- 0.3
Scenario_1_EPS_base$A <- calib_2024$Q_target * (calib_2024$P_target ^ Scenario_1_EPS_base$eps)

path_S1_baseEPS <- simulate_path(
  Scenario_1_EPS_base,
  "Scenario 1",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S1_baseEPS$case <- "Base elasticity"


# --------
# Høy elastisitet
# --------
Scenario_1_EPS_high <- Scenario_1
Scenario_1_EPS_high$eps <- 0.8
Scenario_1_EPS_high$A <- calib_2024$Q_target * (calib_2024$P_target ^ Scenario_1_EPS_high$eps)

path_S1_highEPS <- simulate_path(
  Scenario_1_EPS_high,
  "Scenario 1",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0,
  beta_end = 0,
  alpha_eu_path = rep(0, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)

path_S1_highEPS$case <- "High elasticity"


# SAMLE ELASTISITETSANALYSEN I ÉN TABELL
sens_S1_EPS <- bind_rows(
  path_S1_lowEPS,
  path_S1_baseEPS,
  path_S1_highEPS
)
# ============================================================
# TABELL FOR 2025 OG 2035
# ------------------------------------------------------------
# Henter ut kun startår og sluttår for sensitiviteten
# i Scenario 1 med ulik elastisitet.
# ============================================================

sens_S1_EPS_2025_2035 <- sens_S1_EPS %>%
  filter(year %in% c(2025, 2035)) %>%
  select(case, year, P, Q, x_eu, x_no, x_row) %>%
  arrange(case, year) %>%
  mutate(
    P = round(P, 2),
    Q = round(Q, 2),
    x_eu = round(x_eu, 2),
    x_no = round(x_no, 2),
    x_row = round(x_row, 2)
  )

print(sens_S1_EPS_2025_2035)

# ============================================================
# SENSITIVITETSANALYSE – SCENARIO 2
# ETTERSPØRSELSELASTISITET
# ------------------------------------------------------------
# Vi tester hvordan resultatene i Scenario 2 endres når
# etterspørselselastisiteten varierer mellom lav, base og høy.
# ============================================================

Scenario_2_EPS_low <- Scenario_2A
Scenario_2_EPS_low$eps <- 0.1
Scenario_2_EPS_low$A <- calib_2024$Q_target * (calib_2024$P_target ^ Scenario_2_EPS_low$eps)

path_S2_lowEPS <- simulate_path(
  Scenario_2_EPS_low,
  "Scenario 2",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)
path_S2_lowEPS$case <- "Low elasticity"


Scenario_2_EPS_base <- Scenario_2A
Scenario_2_EPS_base$eps <- 0.3
Scenario_2_EPS_base$A <- calib_2024$Q_target * (calib_2024$P_target ^ Scenario_2_EPS_base$eps)

path_S2_baseEPS <- simulate_path(
  Scenario_2_EPS_base,
  "Scenario 2",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)
path_S2_baseEPS$case <- "Base elasticity"


Scenario_2_EPS_high <- Scenario_2A
Scenario_2_EPS_high$eps <- 0.8
Scenario_2_EPS_high$A <- calib_2024$Q_target * (calib_2024$P_target ^ Scenario_2_EPS_high$eps)

path_S2_highEPS <- simulate_path(
  Scenario_2_EPS_high,
  "Scenario 2",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.69,
  beta_end = 0,
  alpha_eu_path = rep(0.50, length(2025:2035)),
  alpha_no_path = rep(0.42, length(2025:2035)),
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)
path_S2_highEPS$case <- "High elasticity"


sens_S2_EPS <- bind_rows(
  path_S2_lowEPS,
  path_S2_baseEPS,
  path_S2_highEPS
)

sens_S2_EPS_2025_2035 <- sens_S2_EPS %>%
  filter(year %in% c(2025, 2035)) %>%
  select(case, year, P, Q, x_eu, x_no, x_row) %>%
  arrange(case, year) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

print(sens_S2_EPS_2025_2035)









# ============================================================
# LANGSIKTIG SCENARIO 2035
# ------------------------------------------------------------
# Dette scenarioet illustrerer en mulig markedssituasjon i 2035
# der flere sentrale parametere endres samtidig.
#
# Antakelser:
# - høy EU ETS-pris
# - lavere CCS-kostnad
# - høy andel EU-produsenter med CCS
# - moderat subsidieandel
# - ingen gratiskvoter
# - CBAM er aktiv
# ============================================================

Scenario_2035_longrun <- Base_Parameter

# Karbonpris
Scenario_2035_longrun$P_CO2 <- 200

# Gratiskvoter fases helt ut
Scenario_2035_longrun$beta <- 0

# CBAM aktiv
Scenario_2035_longrun$gamma_cbam <- 1
Scenario_2035_longrun$P_home <- 0

# CCS-kostnader faller over tid
Scenario_2035_longrun$C_CCS_eu <- 90
Scenario_2035_longrun$C_CCS_no <- 90

# Subsidieandel
Scenario_2035_longrun$s_eu <- 0.11
Scenario_2035_longrun$s_no <- 0.11

# Andel EU-produsenter med CCS
Scenario_2035_longrun$rho_eu <- 0.78

# Fangstgrader
Scenario_2035_longrun$alpha_eu <- 0.50
Scenario_2035_longrun$alpha_no <- 0.42


# ============================================================
# LANGSIKTIG SCENARIO 2035
# ------------------------------------------------------------
# Dette scenarioet illustrerer en mulig markedssituasjon i 2035
# der flere sentrale parametere endres samtidig.
#
# Antakelser:
# - høy EU ETS-pris
# - lavere CCS-kostnad
# - høy andel EU-produsenter med CCS
# - moderat subsidieandel
# - ingen gratiskvoter
# - CBAM er aktiv
# ============================================================

Scenario_2035_longrun <- Base_Parameter

# Karbonpris
Scenario_2035_longrun$P_CO2 <- 200

# Gratiskvoter fases helt ut
Scenario_2035_longrun$beta <- 0

# CBAM aktiv
Scenario_2035_longrun$gamma_cbam <- 1
Scenario_2035_longrun$P_home <- 0

# CCS-kostnader faller over tid
Scenario_2035_longrun$C_CCS_eu <- 90
Scenario_2035_longrun$C_CCS_no <- 90

# Subsidieandel
Scenario_2035_longrun$s_eu <- 0.11
Scenario_2035_longrun$s_no <- 0.11

# Andel EU-produsenter med CCS
Scenario_2035_longrun$rho_eu <- 0.78

# Fangstgrader
Scenario_2035_longrun$alpha_eu <- 0.50
Scenario_2035_longrun$alpha_no <- 0.42

# ============================================================
# LØS LANGSIKTIG LIKEVEKT I 2035
# ============================================================

eq_2035_longrun <- solve_equilibrium(Scenario_2035_longrun)

cat("\n===== LANGSIKTIG SCENARIO 2035 =====\n")
cat("P:", round(eq_2035_longrun$P, 4), "\n")
cat("Q:", round(eq_2035_longrun$Q, 4), "\n")
cat("x_eu:", round(eq_2035_longrun$x_eu, 4), "\n")
cat("x_no:", round(eq_2035_longrun$x_no, 4), "\n")
cat("x_row:", round(eq_2035_longrun$x_row, 4), "\n")
cat("===================================\n\n")
# ============================================================
# PLOTDATA – SENSITIVITET KARBONPRIS
# ------------------------------------------------------------
# Lager prosentvis endring i P og Q relativt til base case
# for Scenario 1 og Scenario 2 i 2035.
# ============================================================

plot_carbon_sens <- bind_rows(
  sens_S1_carbon %>% mutate(model_scenario = "Scenario 1"),
  sens_S2_carbon %>% mutate(model_scenario = "Scenario 2")
) %>%
  filter(year == 2035) %>%
  filter(case %in% c("Low carbon price", "Base carbon price", "High carbon price")) %>%
  select(model_scenario, case, P, Q) %>%
  pivot_longer(
    cols = c(P, Q),
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(model_scenario, variable) %>%
  mutate(
    base_value = value[case == "Base carbon price"],
    pct_change = 100 * (value - base_value) / base_value
  ) %>%
  ungroup() %>%
  filter(case != "Base carbon price") %>%
  mutate(
    sensitivity_side = case_when(
      case == "High carbon price" ~ "High",
      case == "Low carbon price" ~ "Low"
    ),
    variable = recode(variable,
                      P = "Pris",
                      Q = "Etterspørsel"
    ),
    x_group = paste(sensitivity_side, model_scenario, sep = " - ")
  )

print(plot_carbon_sens)
ggplot(plot_carbon_sens, aes(x = x_group, y = pct_change, fill = variable)) +
  geom_col(
    position = position_dodge(width = 0.65),
    width = 0.5
  ) +
  geom_text(
    aes(
      label = paste0(round(pct_change, 1), "%"),
      y = ifelse(pct_change >= 0, pct_change + 1, pct_change - 1)
    ),
    position = position_dodge(width = 0.65),
    size = 4
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1.5, y = max(plot_carbon_sens$pct_change) + 5, label = "High", size = 5, fontface = "bold") +
  annotate("text", x = 3.5, y = max(plot_carbon_sens$pct_change) + 5, label = "Low", size = 5, fontface = "bold") +
  scale_x_discrete(labels = c(
    "High - Scenario 1" = "S1",
    "High - Scenario 2" = "S2",
    "Low - Scenario 1"  = "S1",
    "Low - Scenario 2"  = "S2"
  )) +
  labs(
    title = "Sensitivitetsanalyse: karbonpris",
    x = NULL,
    y = "Prosentvis endring fra basis (%)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )



# ============================================================
# PLOTDATA – PRODUKSJON (x_eu, x_no, x_row)
# Karbonpris sensitivitet
# ============================================================

plot_prod_sens <- bind_rows(
  sens_S1_carbon %>% mutate(model_scenario = "Scenario 1"),
  sens_S2_carbon %>% mutate(model_scenario = "Scenario 2")
) %>%
  filter(year == 2035) %>%
  filter(case %in% c("Low carbon price", "Base carbon price", "High carbon price")) %>%
  select(model_scenario, case, x_eu, x_no, x_row) %>%
  pivot_longer(
    cols = c(x_eu, x_no, x_row),
    names_to = "region",
    values_to = "value"
  ) %>%
  group_by(model_scenario, region) %>%
  mutate(
    base_value = value[case == "Base carbon price"],
    pct_change = 100 * (value - base_value) / base_value
  ) %>%
  ungroup() %>%
  filter(case != "Base carbon price") %>%
  mutate(
    sensitivity_side = case_when(
      case == "High carbon price" ~ "High",
      case == "Low carbon price" ~ "Low"
    ),
    region = recode(region,
                    x_eu = "EU",
                    x_no = "Norge",
                    x_row = "ROW"
    ),
    x_group = paste(sensitivity_side, model_scenario, sep = " - ")
  )

print(plot_prod_sens)


ggplot(plot_prod_sens, aes(x = x_group, y = pct_change, fill = region)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.55
  ) +
  geom_text(
    aes(
      label = paste0(round(pct_change, 1), "%"),
      y = ifelse(pct_change >= 0, pct_change + 1, pct_change - 1)
    ),
    position = position_dodge(width = 0.7),
    size = 4
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1.5, y = max(plot_prod_sens$pct_change) + 5, label = "High", size = 5, fontface = "bold") +
  annotate("text", x = 3.5, y = max(plot_prod_sens$pct_change) + 5, label = "Low", size = 5, fontface = "bold") +
  scale_x_discrete(labels = c(
    "High - Scenario 1" = "S1",
    "High - Scenario 2" = "S2",
    "Low - Scenario 1"  = "S1",
    "Low - Scenario 2"  = "S2"
  )) +
  labs(
    title = "Sensitivitetsanalyse: karbonpris",
    subtitle = "Prosentvis endring i produksjon fra basis (2035)",
    x = NULL,
    y = "Prosentvis endring (%)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold")
  )




# ============================================================
# PLOTDATA – PRODUKSJON (x_eu, x_no, x_row)
# CCS-kostnad sensitivitet
# ============================================================

plot_prod_ccs <- bind_rows(
  sens_S1_CCS %>% mutate(model_scenario = "Scenario 1"),
  sens_S2_CCS %>% mutate(model_scenario = "Scenario 2")
) %>%
  filter(year == 2035) %>%
  filter(case %in% c("Low CCS cost", "Base CCS cost", "High CCS cost")) %>%
  select(model_scenario, case, x_eu, x_no, x_row) %>%
  pivot_longer(
    cols = c(x_eu, x_no, x_row),
    names_to = "region",
    values_to = "value"
  ) %>%
  group_by(model_scenario, region) %>%
  mutate(
    base_value = value[case == "Base CCS cost"],
    pct_change = 100 * (value - base_value) / base_value
  ) %>%
  ungroup() %>%
  filter(case != "Base CCS cost") %>%
  mutate(
    sensitivity_side = case_when(
      case == "High CCS cost" ~ "High",
      case == "Low CCS cost" ~ "Low"
    ),
    region = recode(region,
                    x_eu = "EU",
                    x_no = "Norge",
                    x_row = "ROW"
    ),
    x_group = paste(sensitivity_side, model_scenario, sep = " - ")
  )

print(plot_prod_ccs)


ggplot(plot_prod_ccs, aes(x = x_group, y = pct_change, fill = region)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.55
  ) +
  geom_text(
    aes(
      label = paste0(round(pct_change, 1), "%"),
      y = ifelse(pct_change >= 0, pct_change + 1, pct_change - 1)
    ),
    position = position_dodge(width = 0.7),
    size = 4
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1.5, y = max(plot_prod_ccs$pct_change) + 5, label = "High", size = 5, fontface = "bold") +
  annotate("text", x = 3.5, y = max(plot_prod_ccs$pct_change) + 5, label = "Low", size = 5, fontface = "bold") +
  scale_x_discrete(labels = c(
    "High - Scenario 1" = "S1",
    "High - Scenario 2" = "S2",
    "Low - Scenario 1"  = "S1",
    "Low - Scenario 2"  = "S2"
  )) +
  labs(
    title = "Sensitivitetsanalyse: CCS-kostnad",
    subtitle = "Prosentvis endring i produksjon fra basis (2035)",
    x = NULL,
    y = "Prosentvis endring (%)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold")
  )


# ============================================================
# PLOTDATA – PRODUKSJON (x_eu, x_no, x_row)
# Elastisitets-sensitivitet
# ============================================================

plot_prod_eps <- bind_rows(
  sens_S1_EPS %>% mutate(model_scenario = "Scenario 1"),
  sens_S2_EPS %>% mutate(model_scenario = "Scenario 2")
) %>%
  filter(year == 2035) %>%
  filter(case %in% c("Low elasticity", "Base elasticity", "High elasticity")) %>%
  select(model_scenario, case, x_eu, x_no, x_row) %>%
  pivot_longer(
    cols = c(x_eu, x_no, x_row),
    names_to = "region",
    values_to = "value"
  ) %>%
  group_by(model_scenario, region) %>%
  mutate(
    base_value = value[case == "Base elasticity"],
    pct_change = 100 * (value - base_value) / base_value
  ) %>%
  ungroup() %>%
  filter(case != "Base elasticity") %>%
  mutate(
    sensitivity_side = case_when(
      case == "High elasticity" ~ "High",
      case == "Low elasticity" ~ "Low"
    ),
    region = recode(region,
                    x_eu = "EU",
                    x_no = "Norge",
                    x_row = "ROW"
    ),
    x_group = paste(sensitivity_side, model_scenario, sep = " - ")
  )

print(plot_prod_eps)


ggplot(plot_prod_eps, aes(x = x_group, y = pct_change, fill = region)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.55
  ) +
  geom_text(
    aes(
      label = paste0(round(pct_change, 1), "%"),
      y = ifelse(pct_change >= 0, pct_change + 1, pct_change - 1)
    ),
    position = position_dodge(width = 0.7),
    size = 4
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1.5, y = max(plot_prod_eps$pct_change) + 5, label = "High", size = 5, fontface = "bold") +
  annotate("text", x = 3.5, y = max(plot_prod_eps$pct_change) + 5, label = "Low", size = 5, fontface = "bold") +
  scale_x_discrete(labels = c(
    "High - Scenario 1" = "S1",
    "High - Scenario 2" = "S2",
    "Low - Scenario 1"  = "S1",
    "Low - Scenario 2"  = "S2"
  )) +
  labs(
    title = "Sensitivitetsanalyse: etterspørselselastisitet",
    subtitle = "Prosentvis endring i produksjon fra basis (2035)",
    x = NULL,
    y = "Prosentvis endring (%)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold")
  )
