// model.js — runs the Random Forest that was trained in Python.
// The notebook (notebook/cub_model.ipynb) exports every scikit-learn tree to model/model-data.js as
// CUB_MODEL.trees. Each tree is a set of arrays indexed by node number:
//   f[node]  feature used to split (-1 means the node is a leaf)
//   t[node]  threshold: go left if feature <= t, otherwise right
//   l / r    index of the left / right child
//   p[node]  for a leaf, the share of training genes in it that were HEG

// Walk one tree from the root to a leaf and return that leaf's P(HEG)
function treeProba(tree, x) {
  let node = 0;
  // scikit-learn compares features as 32-bit floats, so round the same way to get identical answers
  while (tree.f[node] !== -1) node = Math.fround(x[tree.f[node]]) <= tree.t[node] ? tree.l[node] : tree.r[node];
  return tree.p[node];
}

// The forest's answer is the average of all 200 trees (same as scikit-learn's predict_proba)
const forestProba = x => CUB_MODEL.trees.reduce((s, tree) => s + treeProba(tree, x), 0) / CUB_MODEL.trees.length;

// Prediction for one analysed gene: HEG if P(HEG) >= 0.5, otherwise LEG
function predict(r) {
  if (!r.total) return { label: "—", cls: "", text: "Not enough codons", p: NaN };
  const p = forestProba(featureVector(r));
  return p >= 0.5
    ? { label: "HEG", cls: "heg", text: `Highly expressed · P(HEG) = ${p.toFixed(2)}`, p }
    : { label: "LEG", cls: "leg", text: `Lowly expressed · P(HEG) = ${p.toFixed(2)}`, p };
}
