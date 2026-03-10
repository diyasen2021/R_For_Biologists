# R for Life Scientists
### Data Wrangling, Visualisation, and Biological Insight

> *Modern biology generates more data than ever — from field surveys and lab experiments to environmental DNA and genomics. This course teaches you to work confidently with that data in R: how to clean and organise it, explore its structure, and communicate your findings through publication-quality visualisations. No prior coding experience required.*

---

## Who Is This Course For?

This is a **data literacy course for biologists** — not a programming course, and not a bioinformatics course. If you have biological data and need to understand, clean, and visualise it, this course is for you.

It is designed for:

- **Undergraduate and postgraduate biology students** taking their first steps with data analysis beyond Excel
- **Wet lab / molecular biology researchers** working with experimental data who need a more powerful and reproducible workflow
- **Field ecologists and conservation biologists** dealing with species counts, diversity indices, and environmental variables across sites and seasons
- **MSc and PhD researchers** whose supervisors or reviewers expect R-based figures and reproducible analysis scripts

You do **not** need to be a bioinformatician. You do **not** need any prior coding experience. You need curiosity and biological data you want to understand.

---

## What You Will Be Able to Do by the End

- Set up a clean, reproducible R workflow from scratch
- Understand R's core data structures and how biological data maps onto them
- Import, inspect, and clean messy real-world datasets
- Transform and summarise data using `dplyr`
- Build publication-quality figures in `ggplot2` — histograms, density plots, bar charts, box plots, violin plots, and scatter plots
- Apply everything to a self-directed biological data analysis project using a real dataset

---

## Course Contents

The course is delivered as a mixture of **taught lectures and self-paced practicals**. Each lecture is a fully annotated Markdown document with embedded R code that you run directly in RStudio.

| # | Lecture | Key Skills |
|---|---------|------------|
| 01 | **R Syntax and Data Structures** | Comments, variables, `<-`, `=` in functions, vectors, data frames, lists |
| 02 | **Introduction to ggplot2** | Grammar of graphics, layers, aesthetics, geoms, inheritance, themes, `patchwork` |
| — | **Project: Mammalian Sleep, Body Size & Ecology** | End-to-end analysis using the `msleep` dataset |

---

## The Project

The capstone project uses the **`msleep` dataset** — sleep, brain size, body size, and dietary data for 83 mammal species, built directly into `ggplot2`. No file downloads needed. Students apply everything from both lectures to explore real biological hypotheses:

- Do carnivores sleep more than herbivores, and why might that be?
- Does total sleep scale with body size across mammals?
- Does brain size follow the predicted allometric scaling law (brain ∝ body^0.75)?
- Which species have unexpectedly large brains for their body size?

The project has four parts — data cleaning, visualisation, a brain-body allometry challenge (introducing `lm()` and residuals as a stretch component), and an open-ended biological question of the student's own choosing. A full model solution is included.

---

## Requirements

### Software

**1. R** (version 4.0 or higher recommended)
Download from [https://cran.r-project.org](https://cran.r-project.org)

**2. RStudio Desktop** (free version)
Download from [https://posit.co/download/rstudio-desktop](https://posit.co/download/rstudio-desktop)

> **Alternative:** If you prefer not to install anything locally, [Posit Cloud](https://posit.cloud) is a free browser-based version of RStudio that requires no setup.

### R Packages

**For the course (Lectures 01–02 and the Project):**

```r
install.packages(c(
  "ggplot2",           # data visualisation
  "dplyr",             # data wrangling
  "tidyr",             # reshaping data
  "patchwork",         # combining multiple plots
  "ggrepel",           # non-overlapping text labels
  "palmerpenguins"     # penguin dataset used in Lecture 02
))
```

### Prior Knowledge

| You need | You do NOT need |
|---|---|
| Basic familiarity with spreadsheets | Any coding experience |
| An understanding of mean and standard deviation | Statistics beyond introductory level |
| Biological curiosity | A bioinformatics background |

---

## How to Use This Repo

1. Install R and RStudio (see above)
2. Install the required packages for the course
3. Work through the lectures **in order** — later material builds on earlier lectures
4. Open each `.md` file, read the explanatory text, and type the code into a new R Script in RStudio
5. Run each block, observe the output, and experiment — change values and re-run to build intuition
6. Complete the project after finishing both lectures — attempt it before looking at the solution

---

## A Note on Philosophy

This course teaches R as a **biological thinking tool**, not as a programming language to be mastered for its own sake. Every dataset, every plot, and every analysis task has been chosen because it connects to something biologically interesting. The code is the means; the biological insight is the goal.

The best thing you can do in this course is get something wrong, read the error message, fix it, and understand why it works. That is how every working scientist learned to code.

---

## Acknowledgements

- **Palmer Penguins dataset:** Horst AM, Hill AP, Gorman KB (2020). *palmerpenguins: Palmer Archipelago (Antarctica) penguin data.* R package version 0.1.0.
- **msleep dataset:** Savage VM & West GB (2007). A quantitative, theoretical framework for understanding mammalian sleep. *PNAS*, 104(3), 1051–1056.

---

## Licence

All materials in this repository are released under the [MIT Licence](LICENSE) — free to use, adapt, and redistribute with attribution.

---

*Built with R 4.x · ggplot2 · dplyr · tidyverse*
