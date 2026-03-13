library(dplyr)

# Parametere

Base_Parameter <- list(
  ## Etterspørsel
  eps = 0.5,            # ε: priselastisitet i etterspørselen (ε > 0)
  A   = 1e6,            # A: skaleringsparameter for etterspørselen (nivå/markedsstørrelse)
  
  ## Karbonpriser
  P_CO2  = 65,          # P_CO2: karbonpris i EU ETS (€/tCO2)
  P_home = 0,           # P_home: karbonpris i eksportlandet/ROW (€/tCO2)
   
  ## Utslippsintensitet
  I_eu  = 0.80,         # I_eu: utslippsintensitet EU (tCO2 per tonn sement)
  I_no  = 0.80,         # I_no: utslippsintensitet Norge (tCO2 per tonn sement)
  I_row = 0.80,         # I_row: utslippsintensitet ROW (tCO2 per tonn sement)
  
  ## CCS (karbonfangst og -lagring)
  alpha_eu = 0.00,      # α_eu: fangstandel EU (0–1) = andel utslipp fanget med CCS
  alpha_no = 0.0,      # α_no: fangstandel Norge (0–1)
  C_CCS_eu = 84,        # C_CCS_eu: CCS-kostnad EU (€/tonn sement) – fangst+transport+lagring
  C_CCS_no = 84,        # C_CCS_no: CCS-kostnad Norge (€/tonn sement)
  
  ## Gratiskvoter (output-basert tildeling)
  beta = 0.40,  # β: gratistildeling (tCO2 i gratiskvoter per tonn sement) -> verdi: P_CO2*β (€/tonn)
  
 gamma_cbam = 0,   # 0 = ingen CBAM i basisåret, 1 = CBAM på # vet ikke om dette funker men variabel for å 1 ha CBAM på eller 0 ikke CBAM
  
 ## subsidie på CCS
 S_eu = 0.65 * 84,
 S_no = 0.65 * 84,
 
  ## Kostnadsparametere i marginalkostnad: MC_r(x) = C0_r + C1_r*x + (policyledd)
  C0_eu  = 50,          # C0_eu: basekostnad EU (€/tonn) – konstantledd i MC
  C1_eu  = 0.5 ,   # C1_eu: helning EU (€/tonn^2) – økende MC når produksjon øker
  
  C0_no  = 50,          # C0_no: basekostnad Norge (€/tonn)
  C1_no  = 2.5  ,     # C1_no: helning Norge (€/tonn^2)
  
  C0_row = 50,          # C0_row: basekostnad ROW (€/tonn)
  C1_row = 3.0        # C1_row: helning ROW (€/tonn^2)
)

#____________________________________

# Likningene

## Marginal kostnadsfunksjoner MC_x
MC_eu <- function(x, p){ #EU sin MC
  p$C0_eu +
    p$C1_eu * x +
    p$P_CO2 * p$I_eu * (1 - p$alpha_eu) +
    p$alpha_eu * p$C_CCS_eu -
    p$P_CO2 * p$beta -
    p$S_eu * p$alpha_eu
}

MC_no <- function(x, p){ #Norge sin MC
  p$C0_no +
    p$C1_no * x +
    p$P_CO2 * p$I_no * (1 - p$alpha_no) +
    p$alpha_no * p$C_CCS_no -
    p$P_CO2 * p$beta -
    p$S_no * p$alpha_no
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
  P_target     = 160,      # €/tonn
  Q_target     = 327.7,    # Mt
  x_row_target = 11.369,   # Mt
  x_no_target  = 1.132     # Mt
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
  Base_Parameter$alpha_eu * Base_Parameter$C_CCS_eu -
  Base_Parameter$P_CO2 * Base_Parameter$beta -
  Base_Parameter$S_eu * Base_Parameter$alpha_eu


policy_no_2024 <- Base_Parameter$P_CO2 * Base_Parameter$I_no * (1 - Base_Parameter$alpha_no) +
  Base_Parameter$alpha_no * Base_Parameter$C_CCS_no -
  Base_Parameter$P_CO2 * Base_Parameter$beta - 
  Base_Parameter$S_no * Base_Parameter$alpha_no

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
Scenario_REF1$beta <- 0.4
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
Scenario_REF2$beta <- 0.4
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

Scenario_2 <- Base_Parameter

# EU ETS videreføres
Scenario_2$P_CO2 <- 75

# Gratiskvoter settes opp med startverdi og fases ut i simulate_path
Scenario_2$beta <- 0.4

# Startverdi for CCS i EU
# Denne brukes kun som basis, selve banen legges inn under simulate_path
Scenario_2$alpha_eu <- 0

# Norge har CCS
Scenario_2$alpha_no <- 0.42

# Ingen hjemlig karbonpris i ROW
Scenario_2$P_home <- 0

# CBAM styres i simulate_path
Scenario_2$gamma_cbam <- 0


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
# grafer eller videre analyse.

simulate_path <- function(par_base,
                          scenario_name,
                          year_start = 2025,
                          year_end = 2035,
                          beta_start = NULL,
                          beta_end = NULL,
                          alpha_eu_path = NULL,
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
path_REF <- simulate_path(
  Scenario_Reference,
  "Reference",
  year_start = 2025,
  year_end = 2035,
  beta_start = 0.4,
  beta_end = 0,
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
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
# KJØRING AV SCENARIO 2 – EKSEMPEL MED KONSTANT EU-CCS
# ------------------------------------------------------------
# Her antar vi at EU innfører delvis CCS med en konstant gjennomsnittlig
# fangstandel gjennom hele perioden.
#
# Eksempel:
# 30 % av fabrikkene installerer CCS
# og hver av disse fanger 50 % av utslippene
# => gjennomsnittlig alpha_eu = 0.15
# ============================================================

path_S2 <- simulate_path(
  Scenario_2,
  "CBAM+ ccs senario",
  year_start = 2025,
  year_end = 2035,
  
  # gratiskvoter fases ut
  beta_start = 0.4,
  beta_end   = 0,
  
  # eksempel på delvis CCS i EU
  alpha_eu_path = rep(0.15, length(2025:2035)),
  
  # CBAM av i 2025, på fra 2026
  gamma_cbam_path = c(0, rep(1, length(2026:2035)))
)