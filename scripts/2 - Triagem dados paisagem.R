################################################################################
################# Script Network Analysis - Pains Project ######################
################# Dissimilarity in Interaction Networks ########################
################################################################################

### Landscape data screening

# Load tables with landscape data
## Creating colors for plots
col0 = colorRampPalette(c('white', 'cyan', '#007FFF', 'blue', '#00007F'))
col1 = colorRampPalette(c('#7F0000', 'red', '#FF7F00', 'yellow', 'white',
                          'cyan', '#007FFF', 'blue', '#00007F'))
col2 = colorRampPalette(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                          '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                          '#4393C3', '#2166AC', '#053061'))
col3 = colorRampPalette(c('red', 'white', 'blue'))
col4 = colorRampPalette(c('#7F0000', 'red', '#FF7F00', 'yellow', '#7FFF7F',
                          'cyan', '#007FFF', 'blue', '#00007F'))

# Loading land cover data
# 500 m buffer
cobertura_500a <- read.table("solo_5002.csv", sep = ";", header = TRUE, row.names = 1)
cobertura_500a # with correlated variables

# Standardizing the data
cobertura_500a <- decostand(x = cobertura_500a, method = "standardize")

# Correlation
cor500a <- cor(cobertura_500a, method = "pearson")
corrplot::corrplot(cor500a, type = "lower", method = "ellipse", order = "AOE", outline = TRUE, 
                   col = col3(100), cl.length = 21, tl.col = "black", addCoef.col = 'black')

# 1500 m buffer
cobertura_1500 <- read.delim("solo_1500.csv", sep = ";", header = TRUE, row.names = 1)
cobertura_1500 # correlated variables

# Standardizing the data
cobertura_1500 <- decostand(x = cobertura_1500, method = "standardize")

# Correlation
cor1500 <- cor(cobertura_1500, method = "pearson")
corrplot::corrplot(cor1500, type = "lower", method = "ellipse", order = "AOE", outline = TRUE, 
                   col = col3(100), cl.length = 21, tl.col = "black", addCoef.col = 'black')

# 3000 m buffer
cobertura_3000 <- read.delim("solo_30002.csv", sep = ";", header = TRUE, row.names = 1)
cobertura_3000

# Standardizing the data
cobertura_3000 <- decostand(x = cobertura_3000, method = "standardize")

# Correlation
cor3000 <- cor(cobertura_3000, method = "pearson")
corrplot::corrplot(cor3000, type = "lower", method = "ellipse", order = "AOE", outline = TRUE, 
                   col = col3(100), cl.length = 21, tl.col = "black", addCoef.col = 'black')
