setwd("/Users/Michel/Desktop/Research/CRISPR screen/")
theme_set(figure_theme)
setwd("/Users/Michel/Desktop/Research/MFA/code/")

load(file="Akt1_fluxes.Rda")
load(file="D3_1_fluxes.Rda")
load(file="Ras1_fluxes.Rda")

D3_Akt_Ras_fluxes <- rbind(D3_1_fluxes, Akt1_fluxes, Ras1_fluxes)

D3_Akt_Ras_fluxes$data <- factor(D3_Akt_Ras_fluxes$data, levels=c("BMK_D3_1","BMK_Ras_1","BMK_Akt_1"))

D3_Akt_Ras_plot <- ggplot(D3_Akt_Ras_fluxes, aes(x=compound, y=umol.protein.flux,fill=data)) + 
  geom_bar(stat="identity", position="dodge", color="black") + geom_abline(height=0, slope=0) 
D3_Akt_Ras_plot

########
load(file="D3_1_fluxes_no12or20.Rda")
D3_complete_vs_3timepts <- rbind(D3_1_fluxes, D3_1_fluxes_no12or20)
D3_compv3_plot <- ggplot(D3_complete_vs_3timepts, aes(x=compound, y=umol.protein.flux,fill=data)) + 
  geom_bar(stat="identity", position="dodge", color="black") + geom_abline(height=0, slope=0)
D3_compv3_plot

########
load(file="Ras1_fluxes_no24.Rda")
Ras_complete_vs_no24 <- rbind(Ras1_fluxes, Ras1_fluxes_no24)
Ras_comp_no24_plot <- ggplot(Ras_complete_vs_no24, aes(x=compound, y=umol.protein.flux,fill=data)) + 
  geom_bar(stat="identity", position="dodge", color="black") + geom_abline(height=0, slope=0)
Ras_comp_no24_plot