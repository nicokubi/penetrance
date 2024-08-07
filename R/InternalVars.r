#' Default Baseline Data
#'
#' This dataset contains age-specific cancer penetrance rates for both females and males.
#' As an example, the data is derived from the SEER program for colorectal cancer females and males. 
#'
#' @format A data frame with 94 rows and 3 variables:
#' \describe{
#'   \item{Age}{Age in years (1 to 94)}
#'   \item{Female}{Penetrance rate for females}
#'   \item{Male}{Penetrance rate for males}
#' }
baseline_data_default <- data.frame(
  Age = 1:94,
  Female = c(
    0.00000005, 0.00000009, 0.00000021, 0.00000041, 0.00000061, 0.00000081, 0.00000113,
    0.00000225, 0.00000350, 0.00000474, 0.00000599, 0.00000728, 0.00000886, 0.00001049,
    0.00001212, 0.00001375, 0.00001531, 0.00001642, 0.00001746, 0.00001849, 0.00001953,
    0.00002073, 0.00002289, 0.00002521, 0.00002753, 0.00002985, 0.00003262, 0.00003808,
    0.00004400, 0.00004992, 0.00005583, 0.00006215, 0.00007085, 0.00007996, 0.00008906,
    0.00009817, 0.00010812, 0.00012320, 0.00013913, 0.00015506, 0.00017098, 0.00018798,
    0.00021142, 0.00023593, 0.00026044, 0.00028494, 0.00031276, 0.00036047, 0.00041150,
    0.00046251, 0.00051351, 0.00055826, 0.00056563, 0.00056677, 0.00056790, 0.00056903,
    0.00057393, 0.00060153, 0.00063290, 0.00066426, 0.00069559, 0.00072936, 0.00077779,
    0.00082863, 0.00087943, 0.00093019, 0.00098016, 0.00102557, 0.00107018, 0.00111473,
    0.00115923, 0.00120791, 0.00128200, 0.00136022, 0.00143831, 0.00151626, 0.00159429,
    0.00167343, 0.00175258, 0.00183150, 0.00191019, 0.00198761, 0.00205874, 0.00212853,
    0.00219797, 0.00226703, 0.00232237, 0.00229747, 0.00225914, 0.00222074, 0.00218226,
    0.00213527, 0.00203759, 0.00193173
  ),
  Male = c(
    0.00000004, 0.00000009, 0.00000022, 0.00000045, 0.00000068, 0.00000091, 0.00000118,
    0.00000174, 0.00000235, 0.00000296, 0.00000356, 0.00000424, 0.00000539, 0.00000661,
    0.00000783, 0.00000905, 0.00001025, 0.00001134, 0.00001240, 0.00001346, 0.00001452,
    0.00001587, 0.00001894, 0.00002231, 0.00002567, 0.00002904, 0.00003262, 0.00003753,
    0.00004265, 0.00004778, 0.00005290, 0.00005860, 0.00006771, 0.00007738, 0.00008706,
    0.00009673, 0.00010758, 0.00012550, 0.00014460, 0.00016369, 0.00018278, 0.00020326,
    0.00023204, 0.00026221, 0.00029237, 0.00032253, 0.00035783, 0.00042408, 0.00049547,
    0.00056684, 0.00063818, 0.00070320, 0.00073034, 0.00075116, 0.00077195, 0.00079273,
    0.00081720, 0.00086389, 0.00091425, 0.00096456, 0.00101482, 0.00106685, 0.00112977,
    0.00119444, 0.00125904, 0.00132355, 0.00138738, 0.00144751, 0.00150694, 0.00156625,
    0.00162545, 0.00168775, 0.00176926, 0.00185381, 0.00193813, 0.00202223, 0.00210529,
    0.00218325, 0.00226009, 0.00233659, 0.00241273, 0.00248511, 0.00253680, 0.00258475,
    0.00263232, 0.00267951, 0.00271651, 0.00269447, 0.00266248, 0.00263032, 0.00259802,
    0.00255423, 0.00244243, 0.00231964
  )
)

# Unexported internal parameters
#' Currently supported cancer types
CANCER_TYPES <- c(
    "Brain", "Breast", "Colorectal", "Endometrial",
    "Gastric", "Kidney", "Leukemia", "Melanoma", "Ovarian",
    "Osteosarcoma", "Pancreas", "Prostate", "Small Intestine",
    "Soft Tissue Sarcoma", "Thyroid", "Urinary Bladder",
    "Hepatobiliary", "Contralateral"
)

#' Mapping of short and long cancer names
CANCER_NAME_MAP <- list(short = c(
    "BRA", "BC", "COL", "ENDO",
    "GAS", "KID", "LEUK", "MELA", "OC",
    "OST", "PANC", "PROS", "SI",
    "STS", "THY", "UB",
    "HEP", "CBC"
), long = CANCER_TYPES)

#' Currently supported gene types
GENE_TYPES <- c(
    "ATM", "BARD1", "BRCA1", "BRCA2",
    "BRIP1", "CDH1", "CDK4", "CDKN2A", "CHEK2", "EPCAM",
    "MLH1", "MSH2", "MSH6", "MUTYH", "NBN", "PALB2",
    "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53"
)

#' Mapping of genes to supported variants
ALL_GENE_VARIANT_TYPES <- list(
    ATM = "ATM_hetero_anyPV",
    BARD1 = "BARD1_hetero_anyPV",
    BRCA1 = "BRCA1_hetero_anyPV",
    BRCA2 = "BRCA2_hetero_anyPV",
    BRIP1 = "BRIP1_hetero_anyPV",
    CDH1 = "CDH1_hetero_anyPV",
    CDK4 = "CDK4_hetero_anyPV",
    CDKN2A = "CDKN2A[P16]_hetero_anyPV",
    CHEK2 = "CHEK2_hetero_1100delC",
    EPCAM = "EPCAM_hetero_anyPV",
    MLH1 = "MLH1_hetero_anyPV",
    MSH2 = "MSH2_hetero_anyPV",
    MSH6 = "MSH6_hetero_anyPV",
    MUTYH = c(
        "MUTYH_hetero_anyPV",
        "MUTYH_homo_anyPV"
    ),
    NBN = "NBN_hetero_657del5",
    PALB2 = "PALB2_hetero_anyPV",
    PMS2 = "PMS2_hetero_anyPV",
    PTEN = "PTEN_hetero_anyPV",
    RAD51C = "RAD51C_hetero_anyPV",
    RAD51D = "RAD51D_hetero_anyPV",
    STK11 = "STK11_hetero_anyPV",
    TP53 = "TP53_hetero_anyPV"
)