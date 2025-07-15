R package to calculate ecological index of disturbance for pollen assemblage
 ~~~~~~~~~~~~~~~~
 Date (version)  : v0.2.0, 15 Jul 2025
 Author          : Petr Kuneš
 Email           : petr.kunes@natur.cuni.cz
 Citation        : Kuneš, P., Abraham, V., & Herben, T. (2019). 
                   Changing disturbance-diversity relationships in temperate
                   ecosystems over the past 12000 years. Journal of Ecology, 
                   107(4), 1678–1688. doi: 10.1111/1365-2745.13136
 ~~~~~~~~~~~~~~~~


How do you install the package from GitHub?

1. First, you need to install the `devtools` package. Run R and then type


```
install.packages("devtools")
```

2. Load the `devtools` package.

```
library(devtools)
```

3. In most cases, you just use `install_github()` function.


```
install_github("petrkunes/PolEco")
```

4. Load the PolEco package

```
library(PolEco)
```

5. Run the example calculation

```
Eco.Index(PC_Prasilske[1:102,], pollen_plant_dist_2019)
```

