################################################################################
## Project: Desorden Flavivirus
## Script purpose: Convertir datos de salida de mobidb a data.frames
## Date: 1-2-2019 
## Author: FO.
################################################################################

library("ggplot2")
library("reshape2")
library("jsonlite")
library("seqinr")

################################################################################

# La salida de mobidb no es un archivo json valido, hay que editarlo con
# tr '\n' ' ' < IN | sed 's/ } {/} \n{/g' > IN_f.json

#lee el archivo json (editado previamente en bash)
json_file <- "/home/fernando/git/flavivirus-disorder/2019/flavivirus/no_X/hit08/json_in/gp_f.json"
data = stream_in(file(json_file))

#agrega una columna ID incremental
data$id = 1:nrow(data)

#agrega una columna virus con el nombre del virus
lista_virus = data.frame(strsplit(data$acc, "|", fixed = TRUE)) #separa por |
data$virus = gsub("_virus", "", unlist(lista_virus[1,])) #toma la primer fila que son los nombres de los virus

# data frame con id nombre de virus y nombre largo
codigos_virus = data[c("id", "virus", "acc")]

nombres  = data$virus

#subset de trabajo
datos = data[c("id", "virus", "p")]
datos <- datos[-27,] #saco el que estaba mal alineado


salida = "/home/fernando/restaurado/alineamiento/nuevos_scores"

# escribir los datos finales, descomentar

#for (i in 1:nrow(datos)){
#  write(unlist(datos[i, 3]), "/home/fernando/restaurado/alineamiento/nuevos_scores", append=T, sep=",", ncolumns = length(unlist(datos[i, 3])))
#}

################################################################################
# el archivo nuevos_scores es procesado en un script de python para agregar los
# gaps del alineamiento
# /home/fernando/git/flavivirus-disorder/2019/flavivirus/hit08/agrega_gaps_ascores.py
# Con este script se obtiene una matriz de scores con -1 para las posiciones con gaps
################################################################################

# leer los scores con gaps como una matriz 
m = as.matrix(read.table("/home/fernando/restaurado/alineamiento/scores_df2", header = T))

# transformar la matriz de scores en una matriz "binaria" con 1 a las posiciones
# desorneadas y 0 a las ordenadas. (ademas -1 los gaps)
m_binaria = m
m_binaria[!(m_binaria > 0.45)] = "0"
m_binaria[(m_binaria > 0.45)] = "1"
m_binaria[(m == -1)] = "-1"


#para graficar en ggplot es necesario transformar con melt la matriz.

ggplot(melt(m_binaria), aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + scale_fill_manual(values= c("black", "blue4", "red")) + 
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0)) + 
  theme(element_line(size = 0, linetype = 'solid',colour = "black")) + theme(legend.position="none") +
  geom_vline(xintercept = 105, linetype="dotted",  color = "green", size=0.5)


#------------------------------------------------------------------------------#

# Gapstrip ZIKA - Eliminar las columnas donde hay gaps en la secuencia de zika
# fila de zika = 59
s = m[59,] # fila Zika
gaps_zika <- s != -1.00 # gaps en zika
m_gapped <- m[,gaps_zika] # gap strip

#nueva matriz binaria 
m_binaria = m_gapped
m_binaria[!(m_gapped > 0.45)] = "0"
m_binaria[(m_gapped > 0.45)] = "1"
m_binaria[(m_gapped == -1)] = "-1"


# Agregar  lineas verticales (proteinas, zika) lineas horizontales (separador)
zika_lineas <- c(105,215,216,291,791,1143,1369,1499,2116,2243,2266,2517)

df <- data.frame(lineas=zika_lineas)
h_lineas <- 1:59 - 0.5
h_l <- data.frame(hl=h_lineas)
etiquetas = as.character(datos$virus)
etiquetas[10] = "chimeric_sacar"

# graficar 

ggplot(melt(m_binaria), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + scale_fill_manual(values= c("black", "#002089", "red")) + 
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0)) + 
  theme(element_line(size = 0, linetype = 'solid',colour = "black")) + theme(legend.position="none") +
  geom_vline(data=df,aes(xintercept =lineas),linetype="dotted",   color = "gray51", size=0.5) +
  geom_hline(data=h_l,aes(yintercept =h_lineas),  color = "black", size=0.1)+
  scale_y_continuous(breaks= 1:59, labels = etiquetas)

#------------------------------------------------------------------------------#
# sin columnas con gaps 

col_gaps <- c() #guardar columnas con gaps
for(i in 1:nrow(m)){
  gaps <- which(m[i,] == -1)
  col_gaps <- c(col_gaps, gaps)
}
col_gaps <- unique(col_gaps)

m_nogaps <- m[,-col_gaps] # borrar las columnas con gaps de la matriz

#nueva matriz binaria
m_binaria3 =m_nogaps
m_binaria3[!(m_nogaps> 0.45)] = "0"
m_binaria3[(m_nogaps > 0.45)] = "1"

#graficar

ggplot(melt(m_binaria3), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + 
  scale_fill_manual(values= c( "#002089", "red")) + 
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0)) + 
  theme(element_line(size = 0, linetype = 'solid',colour = "black")) + 
  theme(legend.position="none") + scale_y_continuous(breaks= 1:59, labels = etiquetas)

#------------------------------------------------------------------------------#

# grafico con colores continuos sobre la matriz gap stripped zik

# continuos color

jet.colors <- colorRampPalette(c("black", "black", "black", "black", 
                                 "#002089","#002089","#002089","#002089", "#002089", "#002089",
                                 "FCCECE", "red"))

ggplot(melt(m_gapped), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(12)) + 
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0)) + 
  theme(element_line(size = 0, linetype = 'solid',colour = "black")) +   
  geom_vline(data=df,aes(xintercept =lineas),linetype="dotted",   color = "gray51", size=0.5) +
  geom_hline(data=h_l,aes(yintercept =h_lineas),  color = "black", size=0.1)+
  scale_y_continuous(breaks= 1:59, labels = etiquetas)


#------------------------------------------------------------------------------#

# Separar por proteinas (tomando como base la matriz zika gap stripp)

# por ejemplo C, NS3 y NS5
m_c <- m_gapped[,1:105]
m_ns5 <- m_gapped[, 2517:3425]
m_ns3 <- m_gapped[, 1499:2116]

# paleta de colores
jet.colors <- colorRampPalette(c("black", "black", "black", "black", 
                                 "#002089","#002089","#002089","#002089", "#002089", "#002089",
                                 "FCCECE", "red"))
cons.colors <- colorRampPalette(c("yellow", "black", "green"))

ggplot(melt(m_ns3), aes(x=Var2, y=Var1, fill=value)) + geom_tile() + 
  scale_fill_gradientn(colors = jet.colors(12)) + 
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0)) + 
  theme(element_line(size = 0, linetype = 'solid',colour = "black")) +  
  geom_hline(data=h_l,aes(yintercept =h_lineas),  color = "black", size=0.1) +
  scale_y_continuous(breaks= 1:59, labels = etiquetas) 

#------------------------------------------------------------------------------#

# Agregar a los graficos una barra de la conservacion en el alineamiento

archivo_conservacion <- "~/git/flavivirus-disorder/2019/flavivirus/conservacion"

conservacion <- read.table(archivo_conservacion, quote="\"", 
                           comment.char="", 
                           stringsAsFactors=FALSE)

x_ = conservacion$V1
y_ = conservacion$V3


#prueba con una cion que hace un sliding windows para calcular la conservacion
# al parecer puede reemplazarse con un parametro al momento de calcular la C.
k = 20
yo =  rollapply(y_, 2*k-1, function(x) max(rollmean(x, k)), partial = TRUE)

# subset para proteina C
conservacion_C <- conservacion[1:105,]

#graficar matriz completa (zika gap)
ggplot() + geom_tile(data= melt(m_gapped), aes(x=Var2, y=Var1, fill=value)) + 
  scale_fill_gradientn(colors = jet.colors(12)) + 
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0)) + 
  theme(element_line(size = 0, linetype = 'solid',colour = "black")) +   
  geom_hline(data=h_l,aes(yintercept =h_lineas),  color = "black", size=0.1)+
  scale_y_continuous(breaks= 1:59, labels = etiquetas) + 
  geom_tile(data=conservacion, mapping= aes(x=V1, y=63, colour=yo)) + 
  scale_colour_gradientn(colors= cons.colors(3))

#graficar proteina C 
ggplot() + geom_tile(data= melt(m_c), aes(x=Var2, y=Var1)) + 
  scale_fill_gradientn(colors = jet.colors(12)) + 
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0)) + 
  theme(element_line(size = 0, linetype = 'solid',colour = "black")) +   
  geom_hline(data=h_l,aes(yintercept =h_lineas),  color = "black", size=0.1)+
  scale_y_continuous(breaks= 1:59, labels = etiquetas) +  
  geom_tile(data=conservacion_C, mapping= aes(x=V1, y=63, fill=yo)) + 
  scale_fill_gradientn(colors= cons.colors(3))



################################################################################
################################################################################
################################################################################