## Análisis de expresión diferencial con microarreglos en R

### Introducción
El análisis de la expresión de genes mediante la detección de transcriptos específicos, realizado baja una determinada condición experimental, es una de las estrategias más exitosas dentro del campo de la genómica funcional. El conjunto de genes expresados o transcriptos a partir del DNA genómico (transcriptoma) es un determinante fundamental de la función y fenotipo celular. La transcripción es el primer paso hacia la síntesis de proteínas y, por lo tanto, es altamente indicativa de la respuesta celular a los estímulos ambientales o externos. Entender la función de los genes, la dinámica del transcriptoma y el conocimiento del perfil de genes expresados proporcionan una visión sobre los posibles mecanismos regulatorios o bioquímicos, permiten evaluar la posible causa y la consecuencia de las enfermedades, y permiten identificar genes para una intervención farmacológica putativa (Rabinovich, 2004). 
La expresión diferencial es el cambio de los niveles de expresión de uno o más genes entre dos o varios condiciones.

### Objetivos
El objetivo es evaluar el efecto de la variación genética en el cromosoma Y del ratón sobre el tamaño de los cardiomiocitos y la posible dependencia de tales efectos en niveles de testosterona, utilizando datos de perfiles de expresión génica en el tejido cardíaco de ratones. Particularmente, 1) determinar si existe expresión diferencial entre genotipos, 2) determinar si existe expresión diferencial entre tratamientos, y 3) evaluar las diferencias en la respuesta al tratamiento entre los dos genotipos, es decir, si existe una interacción entre los dos genotipos.

### Materiales y métodos

* Diseño experimental
Se ensayaron ocho ratones machos adultos de dos cepas, C57BL/6J y C57BL/6J-chrY<A/J/NaJ>, denominadas B y BY, respectivamente. De cada cepa (genotipo), cuatro animales fueron castrados y cuatro fueron intervenidos con el mismo procedimiento quirúrjuico, excepto que no se realizó la castración (animales intactos usados como control). El ARN se hibridizó a BeadChips Illumina MouseRef-8 v2.0 que contienen ocho microarreglos con 25,697 sondas cada uno. Solo se seleccionaron aleatoriamente 5,000 sondas para este tutorial (ver cómo se realizó más abajo). Para saber más sobre el diseño experimental consultar a Verdugo et al., 2009. 

* Análisis bioinformático
Para realizar el análisis de los datos se descargaron el conjunto de datos completo disponible en la base de datos GEO por ID GSE15354 en http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15354. 
Los scripts y datos necesarios para ejecutar este tutorial se encuentran disponibles en DE_tutorial del servidor genoma.med.uchile.cl de la Universidad de Chile, por lo que fue necesario descargarlos a una computadora personal. 
Los datos se analizaron en el programa RStudio utilizando los comandos que se muestran a continuación:
Descargar e instalar paqueterías disponibles en Bioconductor (http://www.bioconductor.org):
