# Replication of Single-cell Transcriptomics Revealing Inter-ethnic Variation in Immune Response to *Plasmodium falciparum*

This repository contains a **full replication analysis** of the study:

> **Shahin et al. (2025)**  
> *Single-cell transcriptomics reveals inter-ethnic variation in immune response to Falciparum malaria*  
> **The American Journal of Human Genetics**, 112, 709–723

Using publicly available single-cell RNA sequencing data, this project reproduces and visualizes key findings demonstrating how **ethnicity-driven immune programming**, rather than infection status alone, shapes malaria outcomes.

---

## 🎯 Scientific Motivation

Why do Fulani children in West Africa experience **less severe malaria** than neighboring Mossi children despite comparable exposure to *P. falciparum*?

This analysis addresses that question by exploring **baseline and infection-induced immune states** across more than **70,000 single immune cells**, revealing how population-level genetic and transcriptional differences shape disease tolerance.

---

## 🧠 Key Biological Insights (Replication Summary)

### 1️⃣ Ethnicity Dominates Immune Variation  
Across multiple immune cell types, **ethnicity explains 10–40× more transcriptional variation** than infection status itself.

**Implication:**  
Disease severity may be predicted by **pre-infection immune states**, not just pathogen burden.

---

### 2️⃣ Trained Tolerance in Fulani Monocytes  
Fulani monocytes exhibit:
- ↓ Pro-inflammatory cytokines (TNF-α, IL-6, NF-κB)
- ↑ Phagocytic and parasite-clearance receptors (e.g. CD36)

**Interpretation:**  
Parasites are cleared efficiently **without triggering harmful inflammation**, reducing risk of severe malaria.

---

### 3️⃣ Compartmentalized Immune Strategy  
Different immune compartments play distinct roles:

| Cell type | Dominant function |
|---------|------------------|
| Monocytes | Silent parasite clearance |
| B cells | Robust antibody activation |
| CD4⁺ T cells | Dampened activation |
| CD8⁺ / NK cells | Preserved interferon signaling |

**Vaccine insight:**  
Effective vaccines may need to stimulate **B cells without overactivating T cells**.

---

### 4️⃣ Gene-by-Ethnicity Effects  
The same genetic variants behave differently depending on ethnic background:

- **CD36 (rs1049654):** Enhanced phagocytosis specifically in Fulani monocytes  
- **MT2A (rs35402964):** Antioxidant and anti-inflammatory protection

**Conclusion:**  
Genetic effects are **context-dependent**, underscoring the limits of one-size-fits-all precision medicine.

---

### 5️⃣ Metallothioneins as Systemic Protectors  
Widespread upregulation of **MT1X, MT2A** supports:
- Antioxidant defense
- Zinc homeostasis
- Reduced immunopathology
- Enhanced oxidative burst

**Translational angle:**  
Zinc modulation or metallothionein-inducing strategies may be worth exploring.

---

## 🌍 Why This Matters

This work highlights how **systems immunology** and **population genomics** can:
- Explain long-standing clinical observations
- Inform vaccine and therapeutic design
- Address global health inequities
- Elevate under-represented African populations in genomics research

---

## 📊 Dataset & Methods

- **Cells analyzed:** 71,784
- **Participants:** 126 children (Fulani & Mossi)
- **Setting:** Endemic malaria region, Burkina Faso
- **Data source:** GEO (GSE273781, GSE273785)
- **Framework:** Seurat-based scRNA-seq analysis
- **Focus:** Cell-type–specific transcriptional responses to infection and ethnicity

---

## 🧪 Reproducible Analysis Pipeline

Below is the **complete R pipeline** used to generate publication-quality figures replicating the original study.

> **Note:**  
> The pipeline assumes a processed Seurat object (`MoFu_unfiltered.rds`) with metadata fields:
> - `predicted.celltype.l2`
> - `Ethnicity`
> - `Infection`

---
