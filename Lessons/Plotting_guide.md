# Decision Tree  
## “I have X type of data → I should use Y plot”

---

## 1. What is your **main response variable**?

### A. One continuous variable  
*(e.g. abundance, biomass, height, temperature)*

- Histogram / Density plot  
- Boxplot / Violin plot  
- QQ plot *(if planning statistical tests)*

---

### B. One categorical variable  
*(e.g. species, habitat type, treatment)*

- Bar plot *(counts or proportions)*  

---

### C. Presence / Absence data  
*(0/1, detected / not detected)*

- Bar plot of proportions  
- Spatial presence–absence map  
- Time-series presence plot  

---

## 2. Are you **comparing groups**?

### Yes — Continuous response, categorical groups  
*(e.g. abundance across habitats)*

- Boxplots / Violins *(grouped)*  
- Dot / Jitter plots *(small sample sizes)*  

**If more than one grouping variable:**
- Faceted plots  
- Interaction plots  

---

### Yes — Categorical response, categorical groups  

- Stacked bar plots  
- Mosaic plots  

---

### No — Just exploring one group  

- Distribution plots *(histogram, density, boxplot)*  

---

## 3. Are you studying **relationships between variables**?

### Continuous vs Continuous  
*(e.g. temperature vs abundance)*

- Scatter plot  
- Scatter plot with trend line  
- Nonlinear smooth plot  

**If many variables:**
- Pairwise scatterplot matrix  
- Correlation heatmap  

---

### Continuous vs Categorical  
*(e.g. biomass across seasons)*

- Boxplot / Violin plot  
- Faceted histograms  

---

### Categorical vs Categorical  
*(e.g. species × habitat)*

- Stacked bar plot  
- Heatmap of counts or proportions  

---

## 4. Does your data have a **time component**?

### Yes — Single variable over time  

- Line plot  
- Smoothed trend plot  

---

### Yes — Multiple groups over time  

- Grouped line plots  
- Faceted time-series plots  

---

### Yes — Before vs After an event  

- Before–After line plots  
- Difference (delta) plots  

---

## 5. Does your data have a **spatial component**?

### Coordinates (lat/long, grid, transect)

- Spatial scatter plot *(sampling points)*  
- Heatmap / intensity map  
- Faceted spatial maps *(by species or year)*  

---

### Along a gradient  
*(elevation, depth, distance)*

- Transect plots  
- Line plots along gradient  

---

## 6. Are you analysing **communities or species assemblages**?

### Species counts / abundances

- Rank–abundance (Whittaker) plot  
- Species accumulation curve  

---

### Multivariate community structure  

- Ordination plots *(PCA, NMDS, PCoA)*  
- Biplots *(species + environment)*  
- Dendrograms *(clustering)*  

---

## 7. Are you calculating **biodiversity metrics**?

### Single index across groups  

- Bar plots with uncertainty  
- Dot plots with confidence intervals  

---

### Comparing diversity across sampling effort  

- Rarefaction curves  

---

## 8. Are you fitting **models**?

### Model assumptions  

- Residual vs fitted plot  
- QQ plot  
- Scale–location plot  

---

### Model interpretation  

- Effect size plots  
- Coefficient plots  
- Partial dependence plots  

---

### Model validation  

- Observed vs predicted plots  
- Prediction interval plots  

---

## 9. Is your goal **communication or publication**?

### Scientific paper  

- Multi-panel faceted figures  
- Clean boxplots / ordinations  
- Effect-size focused plots  

---

### Policy / conservation / stakeholders  

- Simplified maps  
- Annotated trend plots  
- Clear before–after comparisons  

---

> **Teaching mantra:**  
> *Choose the plot based on the question, not the software.*
