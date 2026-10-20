// app.js — connects the page to features.js and model.js:
// reads the input box, shows the results table, gene details and RSCU chart,
// fills in the model results section and handles the buttons.

const $ = id => document.getElementById(id);
let results = [];
const pct = x => (x * 100).toFixed(1);

// ---------- Details of one gene: summary numbers + RSCU bar chart ----------
function showDetail(r) {
  $("detailName").textContent = r.name;
  $("pred").textContent = r.prediction.label;
  $("predSub").textContent = r.prediction.text;
  $("predCard").className = "stat pred " + r.prediction.cls;
  $("cai").textContent = Number.isNaN(r.cai) ? "—" : r.cai.toFixed(3);
  $("len").textContent = r.total;
  $("lenSub").textContent = r.leftover ? `${r.leftover} trailing base(s) ignored` : "complete codons";
  $("gc").textContent = `${pct(r.gc)} / ${pct(r.gc3)}`;

  // one small group of bars per amino acid; bar length = RSCU (max possible RSCU is 6)
  const chart = $("chart");
  chart.innerHTML = "";
  const maxR = 6;
  const aas = [...new Set(RSCU_CODONS.map(c => CODON_TO_AA[c]))].sort((a, b) => AA_NAMES[a].localeCompare(AA_NAMES[b]));
  for (const aa of aas) {
    const fam = FAMILIES[aa];
    const used = fam.reduce((s, c) => s + r.counts[c], 0);
    const g = document.createElement("div");
    g.className = "aa-group";
    g.innerHTML = `<h5>${AA_NAMES[aa]} <span>(${aa}, n=${used})</span></h5>`;
    for (const c of fam) {
      const v = r.rscu[c];
      const row = document.createElement("div");
      row.className = "bar-row";
      row.title = `${c}: count ${r.counts[c]}, RSCU ${v.toFixed(2)}`;
      row.innerHTML = `<span>${c}</span><div class="bar-track"><div class="bar ${OPTIMAL.has(c) ? "opt" : ""}" style="width:${Math.min(100, (v / maxR) * 100)}%"></div></div><span class="bar-val">${v.toFixed(2)}</span>`;
      g.appendChild(row);
    }
    chart.appendChild(g);
  }
}

// ---------- Results table (shown when more than one sequence is analysed) ----------
function render() {
  $("results").hidden = results.length === 0;
  if (!results.length) return;
  const tbody = $("batch").querySelector("tbody");
  tbody.innerHTML = "";
  $("batchWrap").hidden = results.length < 2;
  results.forEach((r, i) => {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td></td><td>${r.total}</td><td>${pct(r.gc)}</td><td>${pct(r.gc3)}</td><td>${Number.isNaN(r.cai) ? "—" : r.cai.toFixed(3)}</td><td>${Number.isNaN(r.prediction.p) ? "—" : r.prediction.p.toFixed(2)}</td><td class="tag ${r.prediction.cls}">${r.prediction.label}</td>`;
    tr.firstChild.textContent = r.name;          // gene names are set as text, never as HTML
    tr.addEventListener("click", () => {
      tbody.querySelectorAll("tr").forEach(x => x.classList.remove("sel"));
      tr.classList.add("sel");
      showDetail(r);
      $("detail").scrollIntoView({ behavior: "smooth", block: "start" });
    });
    if (i === 0) tr.classList.add("sel");
    tbody.appendChild(tr);
  });
  showDetail(results[0]);
}

// ---------- Analyze button: parse → features → prediction → display ----------
function run() {
  const msg = $("msg");
  msg.className = "msg";
  const records = parseInput($("seq").value);
  if (!records.length) { msg.textContent = "Please paste a DNA sequence first."; results = []; render(); return; }

  const tooShort = records.filter(r => r.seq.length < 3);
  results = records.filter(r => r.seq.length >= 3).map(r => analyze(r.name, r.seq));
  for (const r of results) r.prediction = predict(r);

  const notes = [];
  const invalid = records.reduce((s, r) => s + r.invalid, 0);
  if (invalid) notes.push(`${invalid} non-ACGT character(s) were removed.`);
  if (tooShort.length) notes.push(`${tooShort.length} sequence(s) were too short and skipped.`);
  const brief = results.filter(r => r.total < 100).length;
  if (brief) notes.push(`${brief} sequence(s) have fewer than 100 codons, so the prediction is less reliable.`);
  if (!results.length) { msg.textContent = "No valid sequence found."; render(); return; }

  msg.className = "msg ok";
  msg.textContent = `Analyzed ${results.length} sequence(s). ${notes.join(" ")}`;
  render();
  $("results").scrollIntoView({ behavior: "smooth", block: "start" });
}

// ---------- CSV download of all features ----------
function downloadCsv() {
  const header = ["gene", "codons", "gc", "gc3", "cai", "p_heg", "prediction",
    ...ALL_CODONS.map(c => `freq_${c}`), ...RSCU_CODONS.map(c => `rscu_${c}`)];
  const esc = s => `"${String(s).replace(/"/g, '""')}"`;
  const lines = [header.join(",")];
  for (const r of results) {
    lines.push([esc(r.name), r.total, r.gc.toFixed(4), r.gc3.toFixed(4), Number.isNaN(r.cai) ? "" : r.cai.toFixed(4),
      Number.isNaN(r.prediction.p) ? "" : r.prediction.p.toFixed(4), r.prediction.label,
      ...ALL_CODONS.map(c => r.freq[c].toFixed(5)), ...RSCU_CODONS.map(c => r.rscu[c].toFixed(4))].join(","));
  }
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/csv" }));
  a.download = "codon_features.csv";
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------- Example inputs ----------
// GFP protein written with only E. coli's optimal codons (a synthetic "ideal HEG")
const GFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK";
function backTranslateOptimal(protein) {
  return [...protein].map(aa => FAMILIES[aa].find(c => OPTIMAL.has(c)) || FAMILIES[aa][0]).join("") + "TAA";
}
const wrap60 = s => s.match(/.{1,60}/g).join("\n");

// ---------- Model results section (numbers computed in Python, stored in model/model-data.js) ----------
function renderModel() {
  const m = CUB_MODEL.metrics;
  document.querySelectorAll("[data-k]").forEach(el => { el.textContent = CUB_MODEL[el.dataset.k]; });
  const f = x => (x * 100).toFixed(1) + "%";
  $("mAcc").textContent = f(m.accuracy);
  $("mF1").textContent = m.f1.toFixed(3);
  $("mPrec").textContent = f(m.precision);
  $("mRec").textContent = f(m.recall);
  $("mCai").textContent = f(CUB_MODEL.caiAccuracy);
  $("mCv").textContent = `${f(CUB_MODEL.cvAccuracy[0])} ± ${(CUB_MODEL.cvAccuracy[1] * 100).toFixed(1)}`;
  $("mSplit").textContent = `Trained on ${CUB_MODEL.nTrain} genes, tested on ${CUB_MODEL.nTest} separate genes · ${CUB_MODEL.trees.length} trees · ${CUB_MODEL.trainedWith}`;
  $("cmTP").textContent = m.tp; $("cmFN").textContent = m.fn;
  $("cmFP").textContent = m.fp; $("cmTN").textContent = m.tn;

  const top = CUB_MODEL.importance.map((v, i) => ({ v, ...FEATURE_NAMES[i] }))
    .sort((a, b) => b.v - a.v).slice(0, 15);
  const max = top[0].v;
  $("importance").innerHTML = top.map(t =>
    `<div class="imp-row"><span class="imp-name"><b>${t.codon}</b> ${AA_NAMES[CODON_TO_AA[t.codon]]} · ${t.kind}</span>` +
    `<div class="bar-track"><div class="bar ${OPTIMAL.has(t.codon) ? "opt" : ""}" style="width:${(t.v / max) * 100}%"></div></div>` +
    `<span class="bar-val">${(t.v * 100).toFixed(1)}%</span></div>`).join("");
}

// ---------- Buttons ----------
$("run").addEventListener("click", run);
$("demoGfp").addEventListener("click", () => {
  $("seq").value = ">GFP_optimal_codons\n" + wrap60(backTranslateOptimal(GFP));
  run();
});
$("demoTest").addEventListener("click", () => {
  // held-out NCBI genes the model never saw during training: 4 HEG + 4 LEG
  const pick = label => CUB_MODEL.testGenes.filter(g => g[2] === label).slice(0, 4);
  $("seq").value = [...pick("HEG"), ...pick("LEG")].map(g => `>${g[0]} ${g[1]} (known ${g[2]}: ${g[3]})\n${wrap60(g[4])}`).join("\n");
  run();
});
$("clear").addEventListener("click", () => { $("seq").value = ""; $("msg").textContent = ""; results = []; render(); });
$("csv").addEventListener("click", downloadCsv);
$("file").addEventListener("change", e => {
  const file = e.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => { $("seq").value = reader.result; run(); };
  reader.readAsText(file);
  e.target.value = "";
});

renderModel();
