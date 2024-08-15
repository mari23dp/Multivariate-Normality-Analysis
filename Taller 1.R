#install.packages("tidyverse")
#install.packages("dplyr")
#install.packages("MASS")
#install.packages("knitr")
#install.packages("kableExtra")
#install.packages("car")
require(kableExtra)
require(tidyverse)
require(MASS)
require(MVN)
require(mvShapiroTest)
require(dplyr)
require(knitr)
require(ggplot2)
require(corrplot)
require(psych)
require(car)
require(tinytex)
require(readr)
require(MVA)

datos <- read_delim("anexo1.csv", delim = ";", 
                    escape_double = FALSE, trim_ws = TRUE)

datos <- as.data.frame(lapply(datos, function(x) {
  if(is.character(x)) {
    as.numeric(gsub(",", ".", x))
  } else {
    x
  }
}))


# Verificamos cuántos valores faltantes hay por columna
colSums(is.na(datos))

par(mfrow = c(1, 2))
# Graficamos la matriz de gráficos de dispersión básica
pairs(datos)
par(mfrow = c(1, 1))

##Prueba de normalidad multivariada informal

#Gráfico de dispersión con elipses de confianza

par(mfrow = c(1,3))

# x1 x2
x1x2 <- data.frame(datos$X1, datos$X2)
y1 <- as.matrix.data.frame(x1x2)
bvbox(y1, method = "robust",
      xlab = expression(x[1]),
      ylab = expression(x[2]))

# x1 x3
x1x3 <- data.frame(datos$X1, datos$X3)
y2 <- as.matrix.data.frame(x1x3)
bvbox(y2, method = "robust",
      xlab = expression(x[1]),
      ylab = expression(x[3]))

# x2 x3
x2x3 <- data.frame(datos$X2, datos$X3)
y3 <- as.matrix.data.frame(x2x3)
bvbox(y3, method = "robust",
      xlab = expression(x[2]),
      ylab = expression(x[3]))

# Restablecemos la disposición de gráficos a la configuración predeterminada
par(mfrow = c(1, 1))


### QQplot normal multivariada
X <- as.matrix(datos)
Xbarra<-colMeans(X)
S<-cov(X)
dm <- mahalanobis(X,Xbarra,S)
cuantiles <- qchisq(ppoints(length(X)),df=4)
# Asignamos número de observaciones y dimensión
n <- length(dm)
p <- 3

# Calculamos los cuantiles teóricos de la distribución chi-cuadrado con p grados de libertad
teoricos <- qchisq(ppoints(n), df = p)

# Primero, generamos el gráfico
qqplot(teoricos, dm, main="QQ Plot de las Distancias de Mahalanobis", 
       xlab="Cuantiles Teoricos Chi cuadrado", 
       ylab="Distancias de Mahalanobis Observadas")

# Luego, agregamos la línea diagonal
abline(0, 1, col="red", lwd=2)

# Prueba de normalidad multivariada formal
# ----- Prueba de Shapiro ----- #
require(mvShapiroTest)
mvShapiro.Test(X)
# ----- Otras Pruebas ----- #
mvn(X, mvnTest="mardia") # test de Mardia
mvn(X, mvnTest="hz") # test de Henze-Zirkler
mvn(X, mvnTest="royston") # test de Royston
mvn(X, mvnTest="dh") # test de Doornik-Hansen




#c. Usando la distancia cuadrada generalizada establezca si existen valores atípicos.

par(mfrow = c(1,1.5))

X<-as.matrix(datos)
Xbarra<-colMeans(X)
S<-cov(X)
dm<-mahalanobis(X,Xbarra,S)
cuantiles<-qchisq(ppoints(length(X)),df=4)
Out<-mvn(datos, mvnTest = "mardia",
         multivariateOutlierMethod ="adj") 


# Punto 2
## Realizamos la carga de datos
#install.packages("MultBiplotR")
require(MultBiplotR)

# Renombramos y convertimos a Data Frame
Data <- Protein  %>% as.data.frame()

# Hacemos unas modificaciones para mostrar el nombre del país
Data$Country <- rownames(Data)
rownames(Data) <- NULL  
Data <- Data[, c("Country", colnames(Data)[-ncol(Data)])]

# Mostramos los datos
kable(Data) %>% kable_styling(font_size = 6, full_width = FALSE)


## a. Determine y analice el vector de medias y la matriz de covarianzas muestrales para las **diferentes regiones.**


### Vector De Medias

#Filtramos para quitar Region y Comunist
F_Data <-  Data[, !(colnames(Data) %in% c("Country", "Comunist", "Region"))]

# Calculamos y mostramos el vector de medias
kable(colMeans(F_Data), col.names = "$\\overline{X}$")


### Matriz Covarianzas

# Calculamos y mostramos el vector de medias
kable(cov(F_Data), digits = 3) %>%  kable_styling(font_size = 8, full_width = FALSE)


## b. Calcule la media de las variables por regiones ¿que puede decir al respecto?

# Agrupamos los datos por región
G_Data <- Data %>% group_by(Region)

# Realizamos el calculo de las medias a lo largo de las variables númericas
MediaxRegion <- G_Data %>% summarize(
  across(Red_Meat:Fruits_Vegetables, mean, na.rm = TRUE)
)

#Mostramos las medias por region
kable(MediaxRegion, digits = 3, caption = "Vectores De Medias Por Región") %>%
  kable_styling(font_size = 10, full_width = FALSE)

## c. Intente construir grupos de paises usando representaciones pictóricas (gráficos de estrellas o caras de Chernoff).

# Hacemos la función de escalado de Min-Max
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))}

# Escalamos las variables usando Min-Max 
F_Data_std <- F_Data %>% mutate(across(everything(), min_max_norm)) %>% as.data.frame()    

# Añadimos la columna para indentificar los países
F_Data_std$Country <- Data$Country

# Creamos y mostramos el gráfico de estrellas con ajustes
stars(F_Data_std[, -10], labels = F_Data_std$Country,key.loc = c(15, 1.25), 
      main = "Gráficos de Estrellas para Países", draw.segments = TRUE,
      col.segments = rainbow(ncol(F_Data_std) - 1), cex.main = 1.5, 
      cex.lab = 2, lty = 1.5, ncol = 8)                   

## d. Utilice las herramientas de gráficas más adecuadas para verificar normalidad multivariada.

Data_means<-colMeans(F_Data)
Data_Cov<-cov(F_Data)
dm<-mahalanobis(F_Data ,Data_means,Data_Cov)
cuantiles<-qchisq(ppoints(nrow(F_Data)),df=9)
qqplot(cuantiles,dm)
abline(0,1)


## e. Realice la prueba de Mardia para verificar las hipótesis:
### H0 : Los datos provienen de una poblacion Normal Multivariada
### H1 : Los datos NO provienen de una poblacion Normal Multivariada

# Realizamos los test de normalidad multivariada y univariada
mvn_result <- mvn(F_Data, mvnTest = "mardia", univariateTest = "AD")

# Extraemos los valores que nos interesan para la normalidad multivariada
test_names_mv <- mvn_result$multivariateNormality$Test
p_values_mv <- as.numeric(as.character(mvn_result$multivariateNormality$`p value`))
results_mv <- mvn_result$multivariateNormality$Result

# Extraemos los valores que nos interesan para la normalidad univariada
test_names_uv <- mvn_result$univariateNormality$Test
variables_uv <- mvn_result$univariateNormality$Variable
p_values_uv <- as.numeric(as.character(mvn_result$univariateNormality$`p value`))
results_uv <- mvn_result$univariateNormality$Normality

# Construimos los data frames
results_mv_df <- data.frame(
  Test = test_names_mv,
  Variable = NA,  # Añadimos la columna 'Variable' con NA para mantener la consistencia
  `p-value` = ifelse(is.na(p_values_mv), "NA", round(p_values_mv, 3)),
  Result = results_mv
)

results_uv_df <- data.frame(
  Test = test_names_uv,
  Variable = variables_uv,
  `p-value` = ifelse(is.na(p_values_uv), "NA", round(p_values_uv, 3)),
  Result = results_uv
)

# Renombramos la columna "Result" a "Normality" en ambas tablas para consistencia
colnames(results_mv_df)[4] <- "Result/Normality"
colnames(results_uv_df)[4] <- "Result/Normality"

# Combinamos los data frames para mostrarlos juntos
combined_results_df <- rbind(results_mv_df, results_uv_df)

# Mostramos los resultados en una tabla
kable(combined_results_df, col.names = c("Test", "Variable", "p-value", "Result/Normality")) %>%
  kable_styling(font_size = 12, full_width = FALSE)

## f. Verifique si hay outliers (multivariados) e identifíquelos.

# Realizamos el análisis
Out <- mvn(Data[, c(-1,-2,-3)], mvnTest = "mardia",
           multivariateOutlierMethod = "quan")

# Colocamos los outliers correspondientes con sus países
outlier_labels <- c("Poland", "Greece", "France", "Albania", "Denmark", "Romania", "Portugal", "Spain")
outlier_indices <- c(16, 10, 9, 1, 6, 18, 17, 19) 

# Añadimos la leyenda
legend("bottomright", legend = paste(outlier_indices, outlier_labels, sep = " = "),
       col = "red", pch = 19, cex = 0.8)


### i. De forma **Univariada**

# Obtenemos los nombres de las columnas
variables <- colnames(F_Data) 

# Creamos los valores de mu para las prueba t
mu_value <- c(9, 7, 2, 15, 5, 30, 4, 3, 4)

# Creamos un DataFrame vacío para almacenar los resultados
results <- data.frame(Variable = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Recorremos las columnas y realizamos el t.test para cada una
for (i in 1:length(variables)) {
  variable_name <- variables[i]
  test_result <- t.test(Data[[variable_name]], mu = mu_value[i], conf.level = 0.95)
  
  # Añadimos los resultados al DataFrame
  results <- rbind(results, data.frame(Variable = variable_name, P_Value = test_result$p.value))
}

# Mostramos los resultados de los p-valores
kable(results)

  
  ### ii. De forma **Multivariada**
  
T2_Hot<-function(mu0,alpha,n){
  Data_means<-colMeans(F_Data)
  Data_Cov<-cov(F_Data)
  Inv_Data_Cov <- solve(Data_Cov)
  Dif_Med<- Data_means - mu0
  T_2<-n%*%t(Dif_Med)%*%Inv_Data_Cov%*%Dif_Med
  return(T_2)
}

T2_V <- T2_Hot(mu_value,0.05,25)

alpha=0.05
p=9
n=25
qf<-qf(alpha,p,n-p,lower.tail = F)
Vc<- round((((n-1)*p)/(n-p))*qf,5)

results_df <- data.frame(
  `Valor Estadístico` = T2_V,
  `Valor Crítico` = Vc
)

kable(results_df)

# ----- Prueba de Shapiro ----- #
require(mvShapiroTest)
mvShapiro.Test(X)
# ----- Otras Pruebas ----- #
mvn(X, mvnTest="mardia") # test de Mardia
mvn(X, mvnTest="hz") # test de Henze-Zirkler
mvn(X, mvnTest="royston") # test de Royston
mvn(X, mvnTest="dh") # test de Doornik-Hansen

