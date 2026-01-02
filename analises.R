library(ggplot2)
library(glmmTMB)
library(piecewiseSEM)
library(patchwork)
library(multcomp)
library(ggcorrplot)
library(Hmisc)
library(corrplot)
library(tidyr)
library(dplyr)

data<-read.table("leiothrix_micronutrients.txt", h=T)
summary(data)
data$species<-as.factor(data$species)
data$area<-as.factor(data$area)
data$Amostra<-as.factor(data$Amostra)
data$Planta<-as.factor(data$Planta)

cor.matrix.var<-cor(data[,16:20], method = "pearson")
round(cor.matrix.var,2)

cor.test.var<-cor.mtest(cor.matrix.var, conf.level = 0.95)
cor.test.var


corrplot(cor.matrix.var, p.mat=cor.test.var$p, method='circle', 
         type='lower', insig='blank', addCoef.col='black',
         number.cex=0.8, order = 'AOE', diag = FALSE)

#Micronutrientes não são correlacionados

data_scaled <- data %>%
  mutate(across(c(Zn, Fe, Mn, Cu, B), scale))

#Efeitos nas medidas de Crasifolia----
data_cra<-data_scaled[1:180,]

#Diâmetro----
m0<-glmmTMB(diameter~1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(diameter~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m0,m1,test="F")

summary(m1)
m2<-glmmTMB(diameter~Zn+Fe+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m1,m2,test="F")

summary(m2)
m3<-glmmTMB(diameter~Fe+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m2,m3,test="F")

summary(m3)
#Efeito negativo de Fe e positivo de Cu e B

coefs <- tidy(m3, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Fe", "Cu", "B")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1a<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1a

#Comprimento da folha-----
m0<-glmmTMB(c_folha~1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(c_folha~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m0,m1,test="F")

summary(m1)

m2<-glmmTMB(c_folha~Zn+Fe+Mn+B+ (1|Amostra/area), family=gaussian, data=data_cra)
anova(m1,m2,test="F")

summary(m2)

m3<-glmmTMB(c_folha~Zn+Fe+B+ (1|Amostra/area), family=gaussian, data=data_cra)
anova(m2,m3,test="F")

summary(m3)

m4<-glmmTMB(c_folha~Fe+B+ (1|Amostra/area), family=gaussian, data=data_cra)
anova(m3,m4, test="F")

summary(m4)

#Fe influenciam negativamente e B influencia positivamente o comprimento da folha em L. crassifolia


coefs <- tidy(m4, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Fe", "B")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1b<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1b

#Largura da folha-----
m0<-glmmTMB(l_folha~1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(l_folha~Zn+Fe+Mn+Cu+B + (1|Amostra/area),family=gaussian, data=data_cra)
anova(m0,m1,test="F")

summary(m1)

m2<-glmmTMB(l_folha~Zn+Fe+Mn+B + (1|Amostra/area),family=gaussian, data=data_cra)

anova(m1,m2,test="F")

summary(m2)

m3<-glmmTMB(l_folha~Fe+Mn+B + (1|Amostra/area),family=gaussian, data=data_cra)
anova(m2,m3,test="F")
summary(m3)

#Efeito negativo de Fe e positivo de Mn e B na largura da folha de Crassifolia

coefs <- tidy(m3, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Mn", "Fe", "B")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1c<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1c

#Comprimento da raiz-----

m0<-glmmTMB(c_raiz~1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(c_raiz~Zn+Fe+Mn+Cu+B + (1|Amostra/area), control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")),family=gaussian, data=data_cra)
anova(m0,m1,test="F")

coefs <- tidy(m1, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Zn", "Fe", "Mn", "Cu", "B")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1d<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1d



#Não tem efeito no comprimento da raiz

#Efeitos na biomassa de Crassifollia-----

#Biomassa da Folha----
m0<-glmmTMB(MF~1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(MF~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m0,m1,test="F")

summary(m1)
m2<-glmmTMB(MF~Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m1,m2,test="F")

summary(m2)
m3<-glmmTMB(MF~Fe+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m2,m3,test="F")

summary(m3)

#Efeito negativo de Fe e Positivo de Cu e B ma biomassa de folha

coefs <- tidy(m3, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Cu", "Fe", "B")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1e<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1e


#Biomassa da raiz----
m0<-glmmTMB(MR~ 1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(MR~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m0,m1,test="F")

#Sem efeito de micronutrientes na biomassa de raiz

coefs <- tidy(m1, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Zn", "Cu", "Mn", "Fe", "B")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1f<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1f


#Biomassa do rizoma-----
m0<-glmmTMB(Mriz~ 1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(Mriz~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m0,m1,test="F")

#Sem efeito de micronutrientes na biomassa de rizoma

coefs <- tidy(m1, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Zn", "Cu", "Mn", "Fe", "B")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1g<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1g


#RSRatio----
m0<-glmmTMB(RSRatio~1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(RSRatio~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m0,m1,test="F")

summary(m1)
m2<-glmmTMB(RSRatio~Zn+Fe+Mn+Cu + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m1,m2,test="F")

summary(m2)
m3<-glmmTMB(RSRatio~Zn+Mn+Cu + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m2,m3,test="F")

summary(m3)
m4<-glmmTMB(RSRatio~Zn+Cu + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m3,m4,test="F")

summary(m4)
m5<-glmmTMB(RSRatio~Zn + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m4,m5,test="F")

summary(m5)
#Efeito negativo de Zn na RSRatio

coefs <- tidy(m5, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Zn")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1i<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1i


#SSRatio----
m0<-glmmTMB(SSRatio~1 + (1|Amostra/area), family=gaussian, data=data_cra)
m1<-glmmTMB(SSRatio~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m0,m1,test="F")

summary(m1)
m2<-glmmTMB(SSRatio~Zn+Fe+Mn+Cu + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m1,m2,test="F")

summary(m2)
m3<-glmmTMB(SSRatio~Zn+Fe+Cu + (1|Amostra/area), family=gaussian, data=data_cra)
anova(m2,m3,test="F")

summary(m3)
#Efeito negativo de Zn e Cu e positivo de Fe

coefs <- tidy(m3, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("Zn", "Cu", "Fe")) %>%
  mutate(
    term = factor(term, levels = sort(unique(term)))
  )

fig1j<-ggplot(coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (estimate ± 95% CI)", y = "") +
  theme_classic()

fig1j


#Efeito dos micronutrientes nas medidas de Spiralis-----
data_spi<-data_scaled[181:360,]

#Diâmetro----
m0<-glmmTMB(diameter~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(diameter~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")

#Sem efeito no diâmetro de spiralis
#Comprimento da folha----
m0<-glmmTMB(c_folha~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(c_folha~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")

#Sem efeito no comprimento da folha
#Largura da folha-----
m0<-glmmTMB(l_folha~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(l_folha~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")

summary(m1)

m2<-glmmTMB(l_folha~Zn+Fe+Mn+Cu + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m1,m2,test="F")

summary(m2)

m3<-glmmTMB(l_folha~Zn+Fe+Mn + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m2,m3,test="F")

summary(m3)
#Efeito negativo de Zn e Mn e positivo de Fe

#Efeito no comprimento da raiz-----
m0<-glmmTMB(c_raiz~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(c_raiz~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")

#Sem efeito no comprimento da raiz



#Efeito na biomassa de Spiralis-----
#Biomassa da folha-----
m0<-glmmTMB(MF~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(MF~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")

summary(m1)
m2<-glmmTMB(MF~Zn+Fe+Mn+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m1,m2,test="F")

summary(m2)
m3<-glmmTMB(MF~Zn+Fe+Mn + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m2,m3,test="F")

summary(m3)
#Efeitos negativos de Zn e Mn e positivos de Fe

#Biomassa da raiz-----
m0<-glmmTMB(MR~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(MR~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")

#Sem efeito na biomassa da Raiz
#Biomassa do rizoma-----
m0<-glmmTMB(Mriz~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(Mriz~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")
#Sem efeito na biomassa do rizoma
#RSRatio----
m0<-glmmTMB(RSRatio~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(RSRatio~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")
#Sem efeito na RSRatio
#SSRatio----
m0<-glmmTMB(SSRatio~1 + (1|Amostra/area), family=gaussian, data=data_spi)
m1<-glmmTMB(SSRatio~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_spi)
anova(m0,m1,test="F")
#Sem efeito na SSRatio


#Efeito dos micronutrientes nas medidas de mucronata-----
data_muc<-data_scaled[361:540,]

#Diâmetro----
m0<-glmmTMB(diameter~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(diameter~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")
#Sem efeito no diâmetro de mucronata

#Comprimento da folha-----
m0<-glmmTMB(c_folha~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(c_folha~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")
#Sem efeito no comprimento da folha

#Largura da folha----
m0<-glmmTMB(l_folha~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(l_folha~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")

#Sem efeito na largura da folha
#Comprimento da raiz-----
m0<-glmmTMB(c_raiz~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(c_raiz~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")
#sem efeito no comprimento da raiz


#Efeito na biomassa de mucronata-----
#Biomassa da folha-----
m0<-glmmTMB(MF~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(MF~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")

#Sem efeito na biomassa da folha
#Biomassa da raiz-----
m0<-glmmTMB(MR~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(MR~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")
#Sem efeito na biomassa da raiz
#Biomassa do rizoma-----
m0<-glmmTMB(Mriz~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(Mriz~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")
summary(m1)
m2<-glmmTMB(Mriz~Zn+Fe+Mn+Cu + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m1,m2,test="F")
summary(m2)
m3<-glmmTMB(Mriz~Fe+Mn+Cu + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m2,m3,test="F")
summary(m3)
m4<-glmmTMB(Mriz~Mn+Cu + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m3,m4,test="F")
summary(m4)
#Efeito negativo de Mn e Cu
#RSRatio----
m0<-glmmTMB(RSRatio~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(RSRatio~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")
#Sem efeito na RS Ratio

#SSRatio-----
m0<-glmmTMB(SSRatio~1 + (1|Amostra/area), family=gaussian, data=data_muc)
m1<-glmmTMB(SSRatio~Zn+Fe+Mn+Cu+B + (1|Amostra/area), family=gaussian, data=data_muc)
anova(m0,m1,test="F")
#Sem efeito na SSRatio