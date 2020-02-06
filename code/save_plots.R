## Methods
png(paste0(savewd,"/delay_distribution_backward.png"),height=6,width=8,units="in",res=300)
p_sliding_delays_backward
dev.off()

png(paste0(savewd,"/delay_distribution_forward.png"),height=6,width=8,units="in",res=300)
p_sliding_delays_forward
dev.off()

png(paste0(savewd,"/p_moritz_confirm.png"),height=4,width=7,units="in",res=300)
p_other_confirm_fit
dev.off()

png(paste0(savewd,"/p_incubation.png"),height=4,width=7,units="in",res=300)
p_incubation
dev.off()

png(paste0(savewd,"/p_moritz_hosp.png"),height=4,width=7,units="in",res=300)
p_other_hosp_fit
dev.off()

png(paste0(savewd,"/p_confirm_delay_kudos_global.png"),height=4,width=7,units="in",res=300)
p_confirm_delay_kudos_global
dev.off()

png(paste0(savewd,"/p_confirm_delay_kudos.png"),height=4,width=7,units="in",res=300)
p_confirm_delay_kudos
dev.off()

png(paste0(savewd,"/p_confirm_delay_kudos_gamma.png"),height=4,width=7,units="in",res=300)
p_confirm_delay_kudos_gamma
dev.off()

png(paste0(savewd,"/p_death_delay.png"),height=4,width=7,units="in",res=300)
p_death_delay
dev.off()

## And as pdf
pdf(paste0(savewd,"/delay_distribution_backward.pdf"),height=6,width=8)
p_sliding_delays_backward
dev.off()

pdf(paste0(savewd,"/delay_distribution_forward.pdf"),height=6,width=8)
p_sliding_delays_forward
dev.off()

pdf(paste0(savewd,"/p_moritz_confirm.pdf"),height=4,width=7)
p_other_confirm_fit
dev.off()

pdf(paste0(savewd,"/p_incubation.pdf"),height=4,width=7)
p_incubation
dev.off()

pdf(paste0(savewd,"/p_moritz_hosp.pdf"),height=4,width=7)
p_other_hosp_fit
dev.off()

pdf(paste0(savewd,"/p_confirm_delay_kudos_global.pdf"),height=4,width=7)
p_confirm_delay_kudos_global
dev.off()

pdf(paste0(savewd,"/p_confirm_delay_kudos.pdf"),height=4,width=7)
p_confirm_delay_kudos
dev.off()

pdf(paste0(savewd,"/p_confirm_delay_kudos_gamma.pdf"),height=4,width=7)
p_confirm_delay_kudos_gamma
dev.off()


pdf(paste0(savewd,"/p_death_delay.pdf"),height=4,width=7)
p_death_delay
dev.off()
## Fig 1
png(paste0(savewd,"/main_plot.png"),height=10,width=10,res=300,units="in")
p_result
dev.off()
png(paste0(savewd,"/main_plot_cumu.png"),height=10,width=10,res=300,units="in")
p_result_cumu
dev.off()
pdf(paste0(savewd,"/main_plot.pdf"),height=10,width=10)
p_result
dev.off()
pdf(paste0(savewd,"/main_plot_cumu.pdf"),height=10,width=10)
p_result_cumu
dev.off()

## Fig 2
png(paste0(savewd,"/infections_by_province.png"),height=6,width=12,units="in",res=300)
p_infections
dev.off()
png(paste0(savewd,"/symptoms_by_province.png"),height=6,width=12,units="in",res=300)
p_symptoms
dev.off()

pdf(paste0(savewd,"/infections_by_province.pdf"),height=6,width=12)
p_infections
dev.off()
pdf(paste0(savewd,"/symptoms_by_province.pdf"),height=6,width=12)
p_symptoms
dev.off()

## Fig 3
png(paste0(savewd,"/comparison_plot_geometric.png"),height=6,width=8,units="in",res=300)
plot(p_comparison)
dev.off()

pdf(paste0(savewd,"/comparison_plot_geometric.pdf"),height=6,width=8)
plot(p_comparison)
dev.off()