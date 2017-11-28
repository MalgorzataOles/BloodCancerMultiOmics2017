###########################################################
### ADDITIONAL VARIABLE DEFINITIONS
###########################################################

# Colors for 3 main pathways
pathColor = c("B-cell receptor"="#f12da0",
              "Reactive oxygen species"="#7073b0",
              "Apoptosis (BH3)"="#efa62b",
              "Other"="grey")
              
# Colors for type of drug
typeColor = c("FDA approved"="#23b14d",
              "clinical development/\ntool compound"="#a349a3")

# Assignment of diseases to cell of origin
diagAmt = list(
    `B-cell origin` = c("CLL","MCL","MZL","LPL","B-PLL","HCL","HCL-V","FL"),
    `T-cell origin` = c("T-PLL","Sezary","PTCL-NOS"),
    `Myeloid origin`= c("AML"),
    `Healthy`       = c("hMNC"))
    
# Colors for cell of origin
colDiagL = c("B-cell origin"="skyblue",
    "T-cell origin"="chocolate1",
    "Myeloid origin"="yellowgreen",
    "Healthy"="seashell4")

# Colors for diseases
colDiagS = c(`CLL`="pink", #  WH replaced "violetred" by "pink"
    `MCL`="darkorchid",
    `MZL`="royalblue4",
    `LPL`="cornflowerblue",
    `B-PLL`="darkturquoise",
    `HCL`="cyan4",
    `HCL-V`="mediumaquamarine",
    `FL`="darkseagreen2",
    `T-PLL`="darkgoldenrod1",
    `Sezary`="brown3",
    `PTCL-NOS`="saddlebrown",
    `AML`="olivedrab",
    `hMNC`="black")
    
## colours to annotate sample features (molecular or clinical) in heatmap
sampleColors <- list(
  `del 17p13`  = c(`0`   = "#f8f8f8", `1` = "darkorchid4"),
  del11q22.3   = c(`0`   = "#f8f8f8", `1` = "darkorchid4"),
  del13q14     = c(`0`   = "#f8f8f8", `1` = "darkorchid4"),
  TP53         = c(`0`   = "#f8f8f8", `1` = "darkred"),
  p53          = c(`0`   = "#f8f8f8", `1` = "darkgreen"),
  ATM          = c(`0`   = "#f8f8f8", `1` = "darkgoldenrod3"),
  `trisomy 12` = c(`0`   = "#f8f8f8", `1` = "darkgreen"),
  SF3B1        = c(`0`   = "#f8f8f8", `1` = "darkblue"),
  NOTCH1       = c(`0`   = "#f8f8f8", `1` = "darkgoldenrod3"),
  PRPF8        = c(`0`   = "#f8f8f8", `1` = "darkorchid4"),
  NRAS         = c(`0`   = "#f8f8f8", `1` = "darkgoldenrod3"),
  KRAS         = c(`0`   = "#f8f8f8", `1` = "darkgoldenrod3"),
  `KRAS | NRAS`= c(`0`   = "#f8f8f8", `1` = "darkgoldenrod3"),
  pretreated   = c(`yes` = "firebrick2", `no` = "#f8f8f8"),
  alive        = c(`yes` = "#f8f8f8", `no` = "black"),
  sex          = c(`m`   = "royalblue4", `f` = "maroon"))


# Custom made order of drugs
mainDrOrd = c("D_164","D_030","D_082","D_003","D_002","D_166","D_036","D_079",
              "D_071","D_035","D_050","D_020","D_163","D_078","D_162","D_077",
              "D_075","D_045","D_032","D_074","D_012","D_033","D_023","D_083",
              "D_172","D_029","D_024","D_017","D_063","D_034","D_015","D_021",
              "D_165","D_004","D_040","D_001","D_081","D_048","D_060","D_039",
              "D_067","D_169","D_084","D_025","D_010","D_006","D_159","D_056",
              "D_008","D_127","D_141","D_149","D_049","D_007","D_054","D_013",
              "D_066","D_011","D_168","D_167","D_043","D_041","D_053","D_CHK")
