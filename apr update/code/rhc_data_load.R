# code below borrowed/adapted from: https://github.com/KenLi93/p2sls_surv_manuscript/blob/main/rhc_data_analysis.R
 
wd <- getwd()
 
source(paste0(wd, "/libraries.R"))
 
rhc <- read.csv(paste0(wd,"/data/rhc.csv")) %>%
  mutate(
    ## change the dates to the Date format
    sadmdte = as.Date(sadmdte, origin = "1960-1-1"),
    dschdte = as.Date(dschdte, origin = "1960-1-1"),
    dthdte = as.Date(dthdte, origin = "1960-1-1"),
    lstctdte = as.Date(lstctdte, origin = "1960-1-1"),
    ## dichotomize variables
    death = as.numeric(death == "Yes"),
    sex = as.numeric(sex == "Female"),
    raceblack = as.numeric(race == "black"),
    raceother = as.numeric(race == "other"),
    ## income: under $11k is reference
    income1 = as.numeric(income == "$11-$25k"),
    income2 = as.numeric(income == "$25-$50k"),
    income3 = as.numeric(income == "> $50k"),
    ## insurance type: private is reference
    ins_care = as.numeric(ninsclas == "Medicare"),
    ins_pcare = as.numeric(ninsclas == "Private & Medicare"),
    ins_caid = as.numeric(ninsclas == "Medicaid"),
    ins_no = as.numeric(ninsclas == "No insurance"),
    ins_carecaid = as.numeric(ninsclas == "Medicare & Medicaid"),
    ## primary disease category: ARF is reference
    cat1_copd = as.numeric(cat1 == "COPD"),
    cat1_mosfsep = as.numeric(cat1 == "MOSF w/Sepsis"),
    cat1_mosfmal = as.numeric(cat1 == "MOSF w/Malignancy"),
    cat1_chf = as.numeric(cat1 == "CHF"),
    cat1_coma = as.numeric(cat1 == "Coma"),
    cat1_cirr = as.numeric(cat1 == "Cirrhosis"),
    cat1_lung = as.numeric(cat1 == "Lung Cancer"),
    cat1_colon = as.numeric(cat1 == "Colon Cancer"),
    ## secondary disease category: NA is reference
    cat2_mosfsep = as.numeric(cat2 == "MOSF w/Sepsis" & !is.na(cat2)),
    cat2_coma = as.numeric(cat2 == "Coma" & !is.na(cat2)),
    cat2_mosfmal = as.numeric(cat2 == "MOSF w/Malignancy" & !is.na(cat2)),
    cat2_lung = as.numeric(cat2 == "Lung Cancer" & !is.na(cat2)),
    cat2_cirr = as.numeric(cat2 == "Cirrhosis" & !is.na(cat2)),
    cat2_colon = as.numeric(cat2 == "Colon Cancer" & !is.na(cat2)),
    ## respiratory diagnosis
    resp = as.numeric(resp == "Yes"),
    ## cardiovascular diagnosis
    card = as.numeric(card == "Yes"),
    ## Neurological diagnosis
    neuro = as.numeric(neuro == "Yes"),
    ## Gastrointestinal diagnosis
    gastr = as.numeric(gastr == "Yes"),
    ## Renal diagnosis
    renal = as.numeric(renal == "Yes"),
    ## Metabolic diagnosis
    meta = as.numeric(meta == "Yes"),
    ## Hematologic diagnosis
    hema = as.numeric(hema == "Yes"),
    ## Sepsis diagnosis
    seps = as.numeric(seps == "Yes"),
    ## Trauma diagnosis
    trauma = as.numeric(trauma == "Yes"),
    ## Orthopedic diagnosis
    ortho = as.numeric(ortho == "Yes"),
    ## ADL
    MISSadld3p = as.numeric(is.na(adld3p)),
    ## DNR status on day 1
    dnr1 = as.numeric(dnr1 == "Yes"),
    ## cancer: yes/no/metastatic. No is reference
    ca_yes = as.numeric(ca == "Yes"),
    ca_meta = as.numeric(ca == "Metastatic"),
    # das2d3pc    Duke Activity Status Index, numeric        
    # surv2md1    SUPPORT model estimated survival probability for 2 months
    # aps1        APACHE score, numeric
    # scoma1      Glasgow coma score, numeric
    # wtkilo1     weight in kilo, numeric                                             
    wt0 = as.numeric(wtkilo1 == 0), ## a lot of zero measurements: make extra indicator saying this is 0
    # temp1       temperature, numeric     
    # meanbp1     blood pressure, numeric
    # resp1       respiratory rate, numeric
    # hrt1        heart rate, numeric
    # pafi1       PaO2/FIO2 ratio, numeric
    # paco21      PaCo2, numeric
    # ph1         PH, numeric
    # wblc1       WBC, numeric
    # hema1       Hematocrit, numeric
    # sod1        Sodium, numeric
    # pot1        Potassium, numeric
    # crea1       Creatinine, numeric
    # bili1       Bilirubin, numeric
    # alb1        Albumin, numeric
    # urin1       Urine output, numeric                                                             
    MISSurin1 = as.numeric(is.na(urin1)),  ## ubducatir fir nussubg urine output
    # cardiohx    yes/no, already in 0/1
    # chfhx       yes/no, already in 0/1
    # dementhx    yes/no, already in 0/1
    # psychhx     yes/no, already in 0/1
    # chrpulhx    yes/no, already in 0/1
    # renalhx     yes/no, already in 0/1
    # liverhx     yes/no, already in 0/1
    # gibledhx    yes/no, already in 0/1
    # malighx     yes/no, already in 0/1
    # immunhx     yes/no, already in 0/1
    # transhx     yes/no, already in 0/1
    # amihx       yes/no, already in 0/1
   
    ## Treatment variable
    # swang1      RHC or NO RHC               CATEGORICAL: MAKE BINARY INDICATOR (1)
    swang1 = as.numeric(swang1 == "RHC"),   # 1 indicates RHC was given, 0 means NO RHC
   
    ## time from admission to death
    t2dth = lubridate::time_length(dthdte - sadmdte, unit = "year"),
   
    ## time from last contact to death
    t2lstct = lubridate::time_length(lstctdte - sadmdte, unit = "year"),
    # t2event = ifelse(death == 1, t2dth, t2lstct) + rnorm(length(t2dth)) * 1e-8 # break tie by adding a small random number
    t2event = ifelse(death == 1, t2dth, t2lstct)
  )
