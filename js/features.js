// features.js — turns a DNA sequence into the numbers the model uses.
//   parseInput()  reads plain DNA or FASTA text
//   analyze()     counts codons and computes codon frequency, RSCU, CAI and GC content
//   featureVector() returns the 123 model features (64 frequencies + 59 RSCU values)
// The same calculation is done in Python in notebook/cub_model.ipynb, so the browser
// and the trained model always see identical features.

// ---------- Genetic code ----------
// Codons are listed in the standard T, C, A, G order; AA_STRING gives the amino acid of each one.
const BASES = "TCAG";
const AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const CODON_TO_AA = {};
const ALL_CODONS = [];
{
  let i = 0;
  for (const a of BASES) for (const b of BASES) for (const c of BASES) {
    const codon = a + b + c;
    CODON_TO_AA[codon] = AA_STRING[i++];
    ALL_CODONS.push(codon);
  }
}
const AA_NAMES = {
  A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys", Q: "Gln", E: "Glu", G: "Gly", H: "His", I: "Ile",
  L: "Leu", K: "Lys", M: "Met", F: "Phe", P: "Pro", S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val", "*": "Stop",
};

// Synonymous codon families: amino acid -> the codons that code for it
const FAMILIES = {};
for (const codon of ALL_CODONS) (FAMILIES[CODON_TO_AA[codon]] ||= []).push(codon);

// RSCU is only meaningful for amino acids with more than one codon (no Met, Trp or stops): 59 codons
const RSCU_CODONS = ALL_CODONS.filter(c => CODON_TO_AA[c] !== "*" && FAMILIES[CODON_TO_AA[c]].length > 1);

// ---------- CAI reference: relative adaptiveness w of each codon in E. coli (Sharp & Li 1987) ----------
const W_ECOLI = {
  TTT: 0.296, TTC: 1.000,
  TTA: 0.020, TTG: 0.020, CTT: 0.042, CTC: 0.037, CTA: 0.007, CTG: 1.000,
  ATT: 0.185, ATC: 1.000, ATA: 0.003,
  GTT: 1.000, GTC: 0.066, GTA: 0.495, GTG: 0.221,
  TCT: 1.000, TCC: 0.744, TCA: 0.077, TCG: 0.017, AGT: 0.085, AGC: 0.410,
  CCT: 0.070, CCC: 0.012, CCA: 0.135, CCG: 1.000,
  ACT: 0.965, ACC: 1.000, ACA: 0.076, ACG: 0.099,
  GCT: 1.000, GCC: 0.122, GCA: 0.586, GCG: 0.424,
  TAT: 0.239, TAC: 1.000,
  CAT: 0.291, CAC: 1.000,
  CAA: 0.124, CAG: 1.000,
  AAT: 0.051, AAC: 1.000,
  AAA: 1.000, AAG: 0.253,
  GAT: 0.434, GAC: 1.000,
  GAA: 1.000, GAG: 0.259,
  TGT: 0.500, TGC: 1.000,
  CGT: 1.000, CGC: 0.356, CGA: 0.004, CGG: 0.004, AGA: 0.004, AGG: 0.002,
  GGT: 1.000, GGC: 0.724, GGA: 0.010, GGG: 0.019,
};
// "Optimal" codons (w = 1) are the ones E. coli's highly expressed genes prefer
const OPTIMAL = new Set(Object.keys(W_ECOLI).filter(c => W_ECOLI[c] === 1));

// ---------- Reading the input ----------
// Accepts plain DNA or FASTA (">name" lines). RNA U becomes T; spaces and digits are ignored;
// any other non-ACGT letters are removed and counted so the page can warn about them.
function parseInput(text) {
  text = text.trim();
  if (!text) return [];
  const records = [];
  if (text.startsWith(">")) {
    for (const block of text.split(/^>/m).filter(Boolean)) {
      const lines = block.split(/\r?\n/);
      records.push({ name: lines[0].trim() || `Sequence ${records.length + 1}`, raw: lines.slice(1).join("") });
    }
  } else {
    records.push({ name: "Sequence 1", raw: text });
  }
  return records.map(r => {
    const upper = r.raw.toUpperCase().replace(/U/g, "T").replace(/[\s\d]/g, "");
    const invalid = upper.replace(/[ACGT]/g, "").length;
    return { name: r.name, seq: upper.replace(/[^ACGT]/g, ""), invalid };
  });
}

// ---------- Codon features of one gene ----------
function analyze(name, seq) {
  // 1. count codons, reading in frame from the first base (a trailing partial codon is ignored)
  const counts = Object.fromEntries(ALL_CODONS.map(c => [c, 0]));
  const usable = seq.length - (seq.length % 3);
  for (let i = 0; i < usable; i += 3) counts[seq.slice(i, i + 3)]++;
  const total = usable / 3;

  // 2. codon frequency: share of all codons
  const freq = {};
  for (const c of ALL_CODONS) freq[c] = total ? counts[c] / total : 0;

  // 3. RSCU: observed count ÷ count expected if every synonymous codon were used equally
  const rscu = {};
  for (const fam of Object.values(FAMILIES)) {
    const sum = fam.reduce((s, c) => s + counts[c], 0);
    for (const c of fam) rscu[c] = sum ? counts[c] * fam.length / sum : 0;
  }

  // 4. CAI: geometric mean of the w values of the gene's codons (Met, Trp and stops skipped)
  let logSum = 0, n = 0;
  for (const c of ALL_CODONS) {
    if (!(c in W_ECOLI) || !counts[c]) continue;
    logSum += counts[c] * Math.log(W_ECOLI[c]);
    n += counts[c];
  }
  const cai = n ? Math.exp(logSum / n) : NaN;

  // 5. GC content overall and at the third codon position
  let gc = 0, gc3 = 0;
  for (let i = 0; i < usable; i++) if (seq[i] === "G" || seq[i] === "C") gc++;
  for (let i = 2; i < usable; i += 3) if (seq[i] === "G" || seq[i] === "C") gc3++;

  return {
    name, total, counts, freq, rscu, cai,
    gc: usable ? gc / usable : 0,
    gc3: total ? gc3 / total : 0,
    leftover: seq.length % 3,
  };
}

// The model's 123 inputs, in exactly the order used in the notebook
const FEATURE_NAMES = [...ALL_CODONS.map(c => ({ kind: "Frequency", codon: c })), ...RSCU_CODONS.map(c => ({ kind: "RSCU", codon: c }))];
const featureVector = r => [...ALL_CODONS.map(c => r.freq[c]), ...RSCU_CODONS.map(c => r.rscu[c])];
