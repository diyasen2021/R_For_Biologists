# R Project: Mammalian Sleep, Body Size & Ecology
### An Exploratory Data Analysis in R

---

## Background

You have been given access to a dataset compiled from a landmark study on mammalian sleep patterns published in *Science* (Savage & West, 2007). The dataset — `msleep` — is built into the `ggplot2` package and contains sleep, body size, and dietary data for **83 mammal species**, ranging from the little brown bat to the African elephant.

As a biologist, you know that sleep is not simply a behavioural quirk — it is deeply tied to **metabolic rate**, **predation risk**, **brain development**, and **evolutionary history**. This project asks you to use everything you have learned so far — cleaning, transformation, and visualisation — to explore whether the data support some well-established biological hypotheses.

You will need to think like a biologist, not just a coder. The questions are deliberately open-ended in places, because real data rarely come with neat answers.

---

## The Dataset

Load it with:

```r
library(ggplot2)
library(dplyr)

data(msleep)
head(msleep)
glimpse(msleep)
```

The key variables are:

| Variable | Type | Description |
|---|---|---|
| `name` | Character | Common name of the species |
| `genus` | Character | Taxonomic genus |
| `vore` | Character | Diet: `carnivore`, `herbivore`, `omnivore`, `insectivore` |
| `order` | Character | Taxonomic order (e.g. Primates, Rodentia) |
| `conservation` | Character | IUCN status: `lc`, `nt`, `vu`, `en`, `domesticated` |
| `sleep_total` | Numeric | Total sleep per day (hours) |
| `sleep_rem` | Numeric | REM sleep per day (hours) |
| `sleep_cycle` | Numeric | Length of one sleep cycle (hours) |
| `awake` | Numeric | Hours awake per day |
| `brainwt` | Numeric | Brain weight (kg) |
| `bodywt` | Numeric | Body weight (kg) |

---

## Before You Begin — Explore the Data

Before writing a single line of analysis, always explore your data. Run the following and make a note of what you find.

```r
# How many rows and columns?
dim(msleep)

# What does the structure look like?
str(msleep)

# Summary statistics for every column
summary(msleep)

# How many species per dietary group?
table(msleep$vore)

# How many missing values in each column?
colSums(is.na(msleep))
```

**Question to think about before you start:** Which columns have the most missing data? Why might brain weight be harder to measure than body weight? Does this matter for your analysis?

---

## Part A — Data Cleaning and Transformation
### *Getting the data ready for analysis*

---

### Task A1 — Handle missing values

Not all columns have complete data. For your analyses below, you will need to make deliberate decisions about missing values rather than just dropping everything.

Create two clean subsets of the data:

1. `msleep_clean` — remove rows where `vore` is `NA` (you cannot analyse dietary group without knowing the diet). Keep all other rows, even if other columns have NAs.

2. `msleep_brain` — from `msleep_clean`, keep only rows where **both** `brainwt` and `bodywt` are not `NA`. You will need this for the brain-body allometry task in Part C.

```r
# Your code here
msleep_clean <- ...

msleep_brain <- ...
```

How many species remain in each subset? Use `nrow()` to check.

---

### Task A2 — Create new variables

Body weight in kilograms spans several orders of magnitude in this dataset (from a 0.005 kg shrew to a 6654 kg elephant). When data span this range, a **log transformation** is often essential — it compresses the scale so that patterns across small and large animals are both visible.

Add the following new columns to `msleep_brain` using the `$` operator or `dplyr::mutate()`:

1. `log_bodywt` — the natural log of body weight: `log(bodywt)`
2. `log_brainwt` — the natural log of brain weight: `log(brainwt)`
3. `rem_proportion` — REM sleep as a proportion of total sleep: `sleep_rem / sleep_total`

```r
# Your code here
msleep_brain <- msleep_brain %>%
  mutate(
    log_bodywt     = ...,
    log_brainwt    = ...,
    rem_proportion = ...
  )
```

**Why log-transform?** Look at the raw distribution of `bodywt` with a histogram, then look at the distribution of `log_bodywt`. What changes? Which is more useful for a regression?

---

### Task A3 — Summarise by dietary group

Using `dplyr`, calculate a summary table grouped by `vore`. For each dietary group, calculate:

- the number of species (`n`)
- mean total sleep
- median total sleep
- standard deviation of total sleep
- mean body weight

```r
sleep_summary <- msleep_clean %>%
  group_by(vore) %>%
  summarise(
    n             = n(),
    mean_sleep    = ...,
    median_sleep  = ...,
    sd_sleep      = ...,
    mean_bodywt   = ...,
    .groups = "drop"
  )

sleep_summary
```

Look at the table. Before plotting anything, what pattern do you expect to see? Do carnivores sleep more than herbivores? Why might that be biologically?

---

## Part B — Visualisation
### *Telling the story with plots*

For each plot below, make sure you:

- Add informative axis labels and a title using `labs()`
- Choose a colour palette that works for colour-blind readers (hint: `scale_fill_brewer(palette = "Set2")` or `scale_colour_brewer(palette = "Dark2")` are good choices)
- Use `theme_bw()` or `theme_classic()`

---

### Task B1 — Distribution of total sleep (one numerical variable)

Create a **density plot** showing the distribution of total sleep hours across all species. Colour and fill the curves by `vore` (dietary group). Use `alpha` to make overlapping curves visible.

```r
# Your code here
ggplot(msleep_clean, aes(x = sleep_total, fill = vore, colour = vore)) +
  ...
```

**Interpret your plot:** Is the distribution unimodal or bimodal? Do different dietary groups cluster in different sleep ranges?

---

### Task B2 — Sleep by dietary group (one numerical + one categorical)

Create a **box plot with jittered raw points** showing total sleep broken down by dietary group (`vore`). Colour the boxes by `vore`. Remember to set `outlier.shape = NA` when jitter is present.

```r
# Your code here
```

Now create the same plot but as a **violin plot** with a thin inner box plot. Which version communicates the distribution shape better? Write a one-sentence interpretation beneath each plot as a comment.

---

### Task B3 — Counts per dietary group (one categorical variable)

Create a **bar chart** showing how many species belong to each dietary group. Fill the bars by `vore` and remove the redundant legend.

As an extension, create a second bar chart showing the number of species per `order` (taxonomic order). Use `coord_flip()` to make the long order names readable on the y-axis.

```r
# Extension bar chart
ggplot(msleep_clean, aes(y = order, fill = vore)) +
  geom_bar() +
  ...
```

---

### Task B4 — REM sleep vs total sleep (two numerical variables)

Create a **scatter plot** of `sleep_rem` (y-axis) against `sleep_total` (x-axis). Colour and shape the points by `vore`.

Add a **single overall regression line** (not per group) using `geom_smooth(method = "lm")`. Is there a strong relationship? Does the relationship differ by dietary group?

Then create a second version where you add **separate regression lines per dietary group** by moving `colour = vore` into the global `aes()`. Compare the two plots — does dietary group modify the relationship between REM and total sleep?

```r
# Version 1 — single overall regression line
ggplot(msleep_clean, aes(x = sleep_total, y = sleep_rem)) +
  geom_point(aes(colour = vore, shape = vore), size = 2.5, alpha = 0.8, na.rm = TRUE) +
  ...

# Version 2 — separate lines per group
ggplot(msleep_clean, aes(x = sleep_total, y = sleep_rem, colour = vore)) +
  ...
```

---

### Task B5 — Faceted exploration

Create a **faceted scatter plot** of `log_bodywt` (x-axis) against `sleep_total` (y-axis), faceted by `vore` using `facet_wrap()`. Add a regression line per panel. Add `ncol = 2` to arrange the panels in two columns.

```r
# Your code here
ggplot(msleep_brain, aes(x = log_bodywt, y = sleep_total)) +
  geom_point(aes(colour = vore), size = 2.5, alpha = 0.8) +
  facet_wrap(~ vore, ncol = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  ...
```

**Interpret:** Does the relationship between body size and sleep differ across dietary groups? What does a negative slope in one panel tell you biologically?

---

## Part C — The Brain-Body Allometry Challenge
### *A classic problem in evolutionary biology*

---

### Background: What is allometry?

**Allometry** is the study of how biological traits scale with body size. One of the most famous relationships in biology is that brain size scales predictably with body size across mammals — but not in a simple linear way. The relationship follows a **power law**:

> brain weight ∝ body weight^b

where `b` is the **scaling exponent**. When you log-transform both sides, a power law becomes a straight line:

> log(brain weight) = a + b × log(body weight)

This means that if you plot log brain weight against log body weight, you should see a linear relationship — and the **slope** of that line is the scaling exponent. A value of `b = 1` would mean brain and body scale in direct proportion. In reality for mammals, `b` is around **0.75** — brains scale *sublinearly* with body size (large animals have proportionally smaller brains than small ones).

---

### Task C1 — Plot the allometric relationship

Using `msleep_brain`, create a scatter plot of `log_brainwt` (y-axis) against `log_bodywt` (x-axis). Colour points by `vore`. Add a linear regression line.

```r
# Your code here
ggplot(msleep_brain, aes(x = log_bodywt, y = log_brainwt)) +
  ...
```

Label a few interesting species using `geom_text()` or `ggrepel::geom_text_repel()` (install the `ggrepel` package). For example, label the human, elephant, and little brown bat:

```r
install.packages("ggrepel")
library(ggrepel)

label_species <- c("Human", "African elephant", "Little brown bat", "Horse", "Chimpanzee")

# Filter to just the species you want to label
labels_df <- msleep_brain %>%
  filter(name %in% label_species)

# Add to your plot:
# + geom_text_repel(data = labels_df, aes(label = name), size = 3.5)
```

---

### Task C2 — Fit a linear model

Fit a linear regression of `log_brainwt` on `log_bodywt`. In R, you fit a linear model using `lm()`:

```r
brain_model <- lm(log_brainwt ~ log_bodywt, data = msleep_brain)
summary(brain_model)
```

Read the model output carefully. The `Estimate` for `log_bodywt` is the **slope** — the allometric scaling exponent. 

**Questions to answer:**

1. What is the scaling exponent you estimated? Is it close to the theoretical value of 0.75?
2. What is the R² value (look for `Multiple R-squared`)? What does this tell you about how well body size predicts brain size?
3. Look at the residuals — which species has the largest **positive** residual (bigger brain than expected for its body size)? Which has the largest **negative** residual (smaller brain than expected)? You can extract residuals with `residuals(brain_model)`.

```r
# Find the species with the most extreme residuals
msleep_brain$residual <- residuals(brain_model)

# Largest positive residual (unexpectedly large brain)
msleep_brain %>%
  arrange(desc(residual)) %>%
  select(name, vore, bodywt, brainwt, residual) %>%
  head(5)

# Largest negative residual (unexpectedly small brain)
msleep_brain %>%
  arrange(residual) %>%
  select(name, vore, bodywt, brainwt, residual) %>%
  head(5)
```

Does the species with the largest positive residual surprise you? What might explain an unusually large brain relative to body size?

---

### Task C3 — Does diet group shift the intercept?

Some research suggests that dietary strategy is associated with relative brain size — carnivores and omnivores tend to have larger brains relative to body size than herbivores, possibly because hunting and social complexity demand more neural processing.

Create a new version of the allometry scatter plot, this time adding **separate regression lines per dietary group** using `colour = vore` in the global `aes()`. Do the lines have similar slopes but different intercepts? Or do the slopes differ too?

```r
ggplot(msleep_brain, aes(x = log_bodywt, y = log_brainwt, colour = vore)) +
  geom_point(aes(shape = vore), size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  scale_colour_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(title    = "Brain-Body Allometry by Dietary Group",
       subtitle = "Log-log scale — slope estimates the allometric exponent",
       x        = "log(Body Weight, kg)",
       y        = "log(Brain Weight, kg)",
       colour   = "Diet",
       shape    = "Diet") +
  theme_bw()
```

---

## Part D — Open Exploration
### *Your own biological question*

This is intentionally unstructured. Choose **one** of the following questions (or devise your own), and produce at least **two plots** and a short written interpretation (as comments in your script) to support your answer.

**Option 1 — Predation and sleep:**
The dataset includes a `conservation` column. Animals listed as `lc` (least concern) are often widespread and potentially lower predation risk, while `en` (endangered) or `vu` (vulnerable) species may have very different ecologies. Do animals of different conservation status sleep differently? Is this biological or is it a sampling artefact?

**Option 2 — REM sleep and brain complexity:**
REM sleep is associated with memory consolidation and neural development. Do animals with larger brains (relative to body size — use the residuals you calculated in Task C2) spend a higher proportion of their sleep in REM? Plot `rem_proportion` against your residuals and interpret.

**Option 3 — The sleepiest orders:**
Which taxonomic orders contain the sleepiest mammals? Calculate mean sleep per order and create a ranked horizontal bar chart. Are the top orders what you would expect biologically?

```r
# Option 3 starter
sleep_by_order <- msleep_clean %>%
  group_by(order) %>%
  summarise(mean_sleep = mean(sleep_total, na.rm = TRUE),
            n          = n()) %>%
  filter(n >= 3) %>%           # only orders with at least 3 species
  arrange(desc(mean_sleep))

# Your plot here
ggplot(sleep_by_order, aes(x = mean_sleep, y = reorder(order, mean_sleep))) +
  ...
```

---

## Submission Checklist

Before submitting your project, make sure your R script:

- [ ] Runs from top to bottom without errors on a fresh R session
- [ ] Uses comments (`#`) to explain what each section does and to record your biological interpretations
- [ ] Has meaningful variable names (not `df2`, `plot_final3`, etc.)
- [ ] Uses `<-` for assignment and `=` only inside function arguments
- [ ] Produces all plots from Parts B, C and your chosen Part D question
- [ ] Includes written interpretation of at least three plots as comments

---

## Hints and Further Reading

- **Log transformation:** If you are not sure why we log-transform skewed data, search for "log transformation skewed data R" — there are many good visual explanations.
- **`ggrepel`:** The `ggrepel` package is not covered in lectures but its documentation is very beginner-friendly. Installing and using a new package from CRAN is a skill worth practising.
- **`lm()` output:** The `summary()` of a linear model looks intimidating at first. Focus on the `Estimate` column (the coefficients) and `Multiple R-squared`. The rest will be covered in a later statistics module.
- **Allometry:** For the biological background, see: Jerison, H.J. (1973) *Evolution of the Brain and Intelligence*. Academic Press. A more accessible starting point is the Wikipedia article on *Brain-to-body mass ratio*.
- **`reorder()` in ggplot:** When making ranked bar charts, `reorder(variable, value)` sorts a categorical axis by a numeric variable. This is very useful for ranked comparisons.

---

> *"In biology, nothing makes sense except in the light of evolution — and no pattern becomes visible except in the light of a well-made plot."*

---

## Model Solution

*The solution below is intended for instructors or for self-checking after you have made a genuine attempt. Do not look at it until you have tried each task yourself.*

<details>
<summary><strong>Click to reveal the full solution</strong></summary>

```r
# ============================================================
# PROJECT SOLUTION — Mammalian Sleep, Body Size & Ecology
# ============================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

data(msleep)


# ============================================================
# PART A — Cleaning and Transformation
# ============================================================

# A1 — Remove rows with missing vore
msleep_clean <- msleep %>%
  filter(!is.na(vore))
nrow(msleep_clean)    # 76 species

# Keep only rows with complete brain and body weight
msleep_brain <- msleep_clean %>%
  filter(!is.na(brainwt) & !is.na(bodywt))
nrow(msleep_brain)    # 51 species


# A2 — New variables
msleep_brain <- msleep_brain %>%
  mutate(
    log_bodywt     = log(bodywt),
    log_brainwt    = log(brainwt),
    rem_proportion = sleep_rem / sleep_total
  )


# A3 — Summary table by dietary group
sleep_summary <- msleep_clean %>%
  group_by(vore) %>%
  summarise(
    n            = n(),
    mean_sleep   = round(mean(sleep_total, na.rm = TRUE), 2),
    median_sleep = round(median(sleep_total, na.rm = TRUE), 2),
    sd_sleep     = round(sd(sleep_total, na.rm = TRUE), 2),
    mean_bodywt  = round(mean(bodywt, na.rm = TRUE), 2),
    .groups = "drop"
  )
sleep_summary
# Insectivores and carnivores tend to sleep the most
# Herbivores sleep the least — consistent with higher predation pressure
# and need for longer foraging times


# ============================================================
# PART B — Visualisation
# ============================================================

# B1 — Density plot of total sleep by dietary group
ggplot(msleep_clean, aes(x = sleep_total, fill = vore, colour = vore)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  labs(title   = "Distribution of Total Sleep by Dietary Group",
       x       = "Total Sleep (hours/day)",
       y       = "Density",
       fill    = "Diet",
       colour  = "Diet") +
  theme_bw()
# Insectivores show a notably high-sleep peak
# Herbivores are shifted towards lower sleep hours
# The distribution is somewhat bimodal overall — suggesting two
# broad strategies: short sleepers (large herbivores) and long
# sleepers (small insectivores / carnivores)


# B2 — Box plot with jitter
ggplot(msleep_clean, aes(x = vore, y = sleep_total, fill = vore)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(colour = vore), width = 0.15, size = 2, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  labs(title = "Total Sleep by Dietary Group",
       x     = "Diet",
       y     = "Total Sleep (hours/day)") +
  theme_bw()

# Violin version
ggplot(msleep_clean, aes(x = vore, y = sleep_total, fill = vore)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Total Sleep by Dietary Group — Violin",
       x     = "Diet",
       y     = "Total Sleep (hours/day)") +
  theme_bw()
# The violin reveals a bimodal distribution in herbivores more
# clearly than the box plot — two groups of herbivores with
# different sleep strategies


# B3 — Bar chart: species per dietary group
ggplot(msleep_clean, aes(x = vore, fill = vore)) +
  geom_bar(width = 0.6) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Number of Species per Dietary Group",
       x     = "Diet",
       y     = "Count") +
  theme_bw()

# Extension: species per taxonomic order
ggplot(msleep_clean, aes(y = reorder(order, order, length), fill = vore)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Species per Taxonomic Order",
       y     = "Order",
       x     = "Count",
       fill  = "Diet") +
  theme_bw(base_size = 11)


# B4 — REM vs total sleep
# Version 1: single regression line
ggplot(msleep_clean, aes(x = sleep_total, y = sleep_rem)) +
  geom_point(aes(colour = vore, shape = vore),
             size = 2.5, alpha = 0.8, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "black", linewidth = 0.9, na.rm = TRUE) +
  scale_colour_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(title   = "REM Sleep vs Total Sleep",
       x       = "Total Sleep (hours/day)",
       y       = "REM Sleep (hours/day)",
       colour  = "Diet",
       shape   = "Diet") +
  theme_bw()

# Version 2: separate lines per group
ggplot(msleep_clean,
       aes(x = sleep_total, y = sleep_rem, colour = vore)) +
  geom_point(aes(shape = vore),
             size = 2.5, alpha = 0.8, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE,
              linewidth = 0.9, na.rm = TRUE) +
  scale_colour_brewer(palette = "Set2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(title   = "REM Sleep vs Total Sleep — by Dietary Group",
       x       = "Total Sleep (hours/day)",
       y       = "REM Sleep (hours/day)",
       colour  = "Diet",
       shape   = "Diet") +
  theme_bw()
# Strong positive relationship overall — more total sleep = more REM
# Slopes are broadly similar across groups, suggesting the
# REM:total ratio is relatively conserved across dietary strategies


# B5 — Faceted: body weight vs sleep per dietary group
ggplot(msleep_brain, aes(x = log_bodywt, y = sleep_total)) +
  geom_point(aes(colour = vore), size = 2.5, alpha = 0.8) +
  facet_wrap(~ vore, ncol = 2) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "black", linewidth = 0.7) +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  labs(title    = "Body Size vs Total Sleep by Dietary Group",
       subtitle = "Body weight is log-transformed",
       x        = "log(Body Weight, kg)",
       y        = "Total Sleep (hours/day)") +
  theme_bw()
# Negative relationship in most panels: larger animals sleep less
# Effect is clearest in herbivores and carnivores
# Consistent with the metabolic rate hypothesis: large animals
# have lower mass-specific metabolic rates and need less sleep


# ============================================================
# PART C — Brain-Body Allometry
# ============================================================

# C1 — Allometry scatter plot with species labels
label_species <- c("Human", "African elephant", "Little brown bat",
                   "Horse", "Chimpanzee", "Cow", "Asian elephant")
labels_df <- msleep_brain %>% filter(name %in% label_species)

ggplot(msleep_brain, aes(x = log_bodywt, y = log_brainwt)) +
  geom_point(aes(colour = vore, shape = vore),
             size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey30", linewidth = 0.9) +
  geom_text_repel(data   = labels_df,
                  aes(label = name),
                  size   = 3.2,
                  colour = "grey20") +
  scale_colour_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(title    = "Brain-Body Allometry in Mammals",
       subtitle = "Log-log scale — slope = allometric scaling exponent",
       x        = "log(Body Weight, kg)",
       y        = "log(Brain Weight, kg)",
       colour   = "Diet",
       shape    = "Diet") +
  theme_bw()


# C2 — Linear model
brain_model <- lm(log_brainwt ~ log_bodywt, data = msleep_brain)
summary(brain_model)
# Estimated slope (scaling exponent): ~0.75
# R-squared: ~0.92 — body size explains ~92% of variation in brain size
# This is a remarkably tight relationship across such diverse species

# Examine residuals
msleep_brain$residual <- residuals(brain_model)

# Largest positive residual — bigger brain than expected
msleep_brain %>%
  arrange(desc(residual)) %>%
  select(name, vore, bodywt, brainwt, residual) %>%
  head(5)
# Humans typically have the largest positive residual
# Our brain is ~3x larger than expected for our body size

# Largest negative residual — smaller brain than expected
msleep_brain %>%
  arrange(residual) %>%
  select(name, vore, bodywt, brainwt, residual) %>%
  head(5)
# Large, slow-metabolising animals often have relatively small brains


# C3 — Allometry by dietary group
ggplot(msleep_brain,
       aes(x = log_bodywt, y = log_brainwt, colour = vore)) +
  geom_point(aes(shape = vore), size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  scale_colour_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(title    = "Brain-Body Allometry by Dietary Group",
       subtitle = "Parallel slopes suggest similar scaling; offset = different brain investment",
       x        = "log(Body Weight, kg)",
       y        = "log(Brain Weight, kg)",
       colour   = "Diet",
       shape    = "Diet") +
  theme_bw()
# Slopes are broadly similar across dietary groups (~0.75)
# Intercepts differ: carnivores and omnivores sit above herbivores
# This means at the same body size, carnivores have larger brains
# Consistent with the cognitive demands of hunting and social living


# ============================================================
# PART D — Open Exploration: Option 3
# Sleepiest taxonomic orders
# ============================================================

sleep_by_order <- msleep_clean %>%
  group_by(order) %>%
  summarise(
    mean_sleep = mean(sleep_total, na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  ) %>%
  filter(n >= 2) %>%
  arrange(desc(mean_sleep))

ggplot(sleep_by_order,
       aes(x = mean_sleep,
           y = reorder(order, mean_sleep),
           fill = mean_sleep)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("n=", n)),
            hjust = -0.2, size = 3, colour = "grey40") +
  scale_fill_gradient(low = "#A8D5E2", high = "#1F6B8E",
                      guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Mean Total Sleep by Taxonomic Order",
       subtitle = "Orders with ≥ 2 species | n = number of species",
       x        = "Mean Total Sleep (hours/day)",
       y        = "Taxonomic Order") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank())
# Chiroptera (bats) are the sleepiest order — consistent with
# their high metabolic cost of flight and echolocation
# Proboscidea (elephants) sleep the least — large herbivores
# with low mass-specific metabolic rate and high foraging demands
# Primates are in the middle — moderate sleep with high REM
```

</details>
