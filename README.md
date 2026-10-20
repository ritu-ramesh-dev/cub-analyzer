# CUB Analyzer

**Machine Learning-Based Analysis of Codon Usage Bias to Distinguish Highly and Lowly Expressed Genes in Prokaryotes**

IBD Project · Group A09 — Ritu Ramesh, Dinesh Karthick Raja, Jai Akash, Jayapriya

A working model website: paste an *E. coli* coding sequence (or upload a FASTA file) and the page calculates codon frequency, RSCU, GC / GC3 and CAI, and a Random Forest trained in Python predicts whether the gene is highly expressed (HEG) or lowly expressed (LEG).

## Running it

- **Website:** open `index.html` in any browser (double-click it). It needs no server and no internet connection.
- **Model:** `pip install -r notebook/requirements.txt`, then run `notebook/cub_model.ipynb` in Jupyter (or Google Colab, after `!pip install biopython`). Running it re-creates `model/model-data.js` and `data/cub_dataset.csv`.

## How the code works

```
notebook/cub_model.ipynb  ──trains──►  model/model-data.js  ──loaded by──►  index.html
   (Python)                              (the 200 trees)                      js/features.js → js/model.js → js/app.js
```

| File | What it does |
|---|---|
| `notebook/cub_model.ipynb` | Python pipeline: downloads the genes from NCBI (Biopython), labels HEG / LEG, builds the 123 codon features (pandas), splits train / test, trains and evaluates the Random Forest (scikit-learn), plots the results (matplotlib / seaborn) and exports the model |
| `model/model-data.js` | The trained forest exported by the notebook: every tree's split features, thresholds and leaf probabilities, plus the test metrics, feature importances and the held-out test genes |
| `js/features.js` | Reads plain DNA or FASTA, counts codons and computes codon frequency, RSCU, CAI and GC content — the same calculation as in the notebook |
| `js/model.js` | Runs the exported trees: each tree is walked from root to leaf, and the forest's P(HEG) is the average of the 200 leaves |
| `js/app.js` | Page logic: Analyze button, results table, per-gene RSCU chart, CSV download, examples, and the model results section |
| `index.html` | Page structure: analyzer, how it works, model results, team |
| `css/style.css` | Styling |
| `data/cub_dataset.csv` | The 208 labelled genes with all 123 features and their train / test split |
| `data/ecoli_k12_cds.fna` | The NCBI download (all 4,318 *E. coli* K-12 coding sequences) |

## Data and method

1. **Data:** all 4,318 coding sequences of *E. coli* K-12 MG1655 from NCBI RefSeq NC_000913.3; duplicate records of the same gene are removed.
2. **Labels** from NCBI gene annotation:
   - **104 HEG:** ribosomal proteins (except the rarely made paralogs ykgM / ykgO), translation factors, chaperones, major outer-membrane proteins, glycolysis, TCA and ATP-synthase genes
   - **104 LEG:** DNA-binding transcriptional regulators (highly expressed global regulators excluded), evenly sampled from 252 candidates
3. **Features:** 64 codon frequencies + 59 RSCU values per gene = 123 features.
4. **Train / test split:** `StratifiedGroupKFold` holds out 25% of the genes (156 train / 52 test, balanced). Near-identical gene copies (e.g. *tufA* / *tufB*) are kept on the same side, and the notebook checks that no gene or overlapping sequence is shared.
5. **Model:** `RandomForestClassifier` with 200 trees, trained on the training genes only.

## Results (52 held-out test genes)

| Metric | Value |
|---|---|
| Accuracy | 98.1% |
| F1-score | 0.980 |
| Precision / Recall (HEG) | 100% / 96.2% |
| Grouped 5-fold cross-validation accuracy | 95.7% ± 2.8% |
| CAI rule (CAI ≥ 0.50 = HEG), same test genes | 98.1% |

Cross-validation averages five different train / test splits, so it is the more reliable estimate. The most important codons (CTG, CGT, AAC, GGT, ...) are the classic preferred codons of highly expressed *E. coli* genes.
