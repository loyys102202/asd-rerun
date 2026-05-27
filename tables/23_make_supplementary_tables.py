"""Same three-line treatment for Supplementary Tables."""
import json, pandas as pd, numpy as np
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

def remove_all_borders(table):
    tbl = table._tbl
    tblPr = tbl.find(qn('w:tblPr'))
    if tblPr is None:
        tblPr = OxmlElement('w:tblPr'); tbl.insert(0, tblPr)
    existing = tblPr.find(qn('w:tblBorders'))
    if existing is not None: tblPr.remove(existing)
    tblBorders = OxmlElement('w:tblBorders')
    for edge in ['top','left','bottom','right','insideH','insideV']:
        b = OxmlElement(f'w:{edge}'); b.set(qn('w:val'), 'nil'); tblBorders.append(b)
    tblPr.append(tblBorders)
    for row in table.rows:
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            ex = tcPr.find(qn('w:tcBorders'))
            if ex is not None: tcPr.remove(ex)
            tcB = OxmlElement('w:tcBorders')
            for edge in ['top','left','bottom','right']:
                b = OxmlElement(f'w:{edge}'); b.set(qn('w:val'), 'nil'); tcB.append(b)
            tcPr.append(tcB)

def set_row_border(table, row_idx, edge='bottom', sz=12):
    for cell in table.rows[row_idx].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tcBorders = tcPr.find(qn('w:tcBorders'))
        if tcBorders is None:
            tcBorders = OxmlElement('w:tcBorders'); tcPr.append(tcBorders)
        existing = tcBorders.find(qn(f'w:{edge}'))
        if existing is not None: tcBorders.remove(existing)
        b = OxmlElement(f'w:{edge}')
        b.set(qn('w:val'), 'single'); b.set(qn('w:sz'), str(sz)); b.set(qn('w:color'), '000000')
        tcBorders.append(b)

def remove_cell_shading(cell):
    tcPr = cell._tc.get_or_add_tcPr()
    ex = tcPr.find(qn('w:shd'))
    if ex is not None: tcPr.remove(ex)
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear'); shd.set(qn('w:color'), 'auto'); shd.set(qn('w:fill'), 'auto')
    tcPr.append(shd)

def make_three_line_table(doc, rows_data, font_size=9):
    n_rows = len(rows_data); n_cols = len(rows_data[0])
    t = doc.add_table(rows=n_rows, cols=n_cols)
    t.style = 'Table Grid'; t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for ri, row_data in enumerate(rows_data):
        for ci, txt in enumerate(row_data):
            cell = t.cell(ri, ci)
            cell.text = ""
            p = cell.paragraphs[0]
            run = p.add_run(str(txt))
            run.font.size = Pt(font_size); run.font.name = 'Times New Roman'
            if ri == 0: run.font.bold = True
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            remove_cell_shading(cell)
    remove_all_borders(t)
    set_row_border(t, 0, 'top', 12)
    set_row_border(t, 0, 'bottom', 6)
    set_row_border(t, n_rows-1, 'bottom', 12)
    return t

def add_caption(doc, text, size=10, space_before=Pt(14)):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = space_before; p.paragraph_format.space_after = Pt(4)
    run = p.add_run(text); run.font.size = Pt(size); run.font.bold = True
    run.font.name = 'Times New Roman'
    return p

# Load artifacts
with open('../locked_model_main_8gene.json') as f: m8 = json.load(f)
with open('../locked_model_sens_14gene.json') as f: m14 = json.load(f)
auc_by_rep = pd.read_csv('../nestedcv_auc_by_rep.csv')
n_genes = pd.read_csv('../nestedcv_n_genes.csv')
with open('../cell_marker_analysis.json') as f: cm = json.load(f)
with open('../cell_marker_GSE42133.json') as f: cm42 = json.load(f)
with open('../GSE25507_exploration.json') as f: g25 = json.load(f)
ph = {key: pd.read_csv(f'../phenotype_{key}.csv') for key in ['GPL570','GPL6244','GSE42133','GSE25507']}

# ---- Doc ----
doc = Document()
section = doc.sections[0]
section.top_margin = Cm(2.0); section.bottom_margin = Cm(2.0)
section.left_margin = Cm(1.8); section.right_margin = Cm(1.8)
style = doc.styles['Normal']; style.font.name = 'Times New Roman'; style.font.size = Pt(10)

title = doc.add_paragraph(); title.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = title.add_run("Supplementary Tables — Revised manuscript ISCIENCE-D-26-03472")
r.font.size = Pt(13); r.font.bold = True; r.font.name = 'Times New Roman'

# Table S1
add_caption(doc, "Table S1.  Sensitivity analyses for the LASSO + stability-"
                  "selection panel at alternative thresholds and inclusion "
                  "rules.  The 8-gene panel collapses to 7 unique genes after "
                  "dropping the perfectly co-linear CES1P1 pseudogene.  The "
                  "14-gene panel at a relaxed 50 % stability threshold produces "
                  "multicollinear inflated coefficients and an unrealistically "
                  "high apparent AUC of 1.000; it is reported only to "
                  "demonstrate the failure mode that the 60 % threshold was "
                  "chosen to avoid.", space_before=Pt(8))
s1 = [["Gene","Biotype","8-gene panel coefficient","14-gene panel coefficient","Notes"]]
all_genes = list(dict.fromkeys(m8['genes'] + m14['genes']))
m8_coef = dict(zip(m8['genes'], m8['coefficients']))
m14_coef = dict(zip(m14['genes'], m14['coefficients']))
biotype_map = {
    'KYNU':'Protein-coding','CMYA5':'Protein-coding','MAPK8IP1':'Protein-coding',
    'CES1':'Protein-coding','CES1P1':'Pseudogene (CES1)',
    'KIAA0087':'lncRNA','POU2AF1':'Protein-coding (B-cell)',
    'ERVK3-2':'Endogenous retroviral','CALCOCO1':'Protein-coding',
    'PCDHB2':'Protein-coding','GABRB1':'Protein-coding',
    'LOC728392':'Uncharacterised','STXBP6':'Protein-coding','CDC42':'Protein-coding',
}
notes_map = {
    'CES1P1':'Co-linear with CES1; same coefficient',
    'CDC42':'Adds noise at 50 % threshold',
    'STXBP6':'Adds noise at 50 % threshold',
    'CALCOCO1':'Only stable at 50 % threshold',
    'PCDHB2':'Only stable at 50 % threshold',
    'GABRB1':'Only stable at 50 % threshold',
    'LOC728392':'Only stable at 50 % threshold',
}
for g in all_genes:
    c8 = f"{m8_coef[g]:+.4f}" if g in m8_coef else "—"
    c14 = f"{m14_coef[g]:+.4f}" if g in m14_coef else "—"
    s1.append([g, biotype_map.get(g,''), c8, c14, notes_map.get(g,'')])
s1.append(["Intercept (β₀)","—",f"{m8['intercept']:+.4f}",f"{m14['intercept']:+.4f}","—"])
s1.append(["Apparent AUC on training","—",f"{m8['apparent_auc']:.4f}",
            f"{m14['apparent_auc']:.4f}","Compare to locked 7-gene 0.9421"])
s1.append(["Max |pairwise correlation| of selected genes","—",
            f"{m8['max_abs_corr']:.4f}",f"{m14['max_abs_corr']:.4f}",
            "1.000 → multicollinear"])
make_three_line_table(doc, s1, font_size=8)

# Table S2
add_caption(doc, "Table S2.  Distribution of per-repeat mean AUCs across 100 "
                  "repeats of 5-fold stratified nested cross-validation, with "
                  "outer-fold and per-fold gene-count summary statistics.")
q025 = np.percentile(auc_by_rep['auc_rep_mean'], 2.5)
q975 = np.percentile(auc_by_rep['auc_rep_mean'], 97.5)
mean = auc_by_rep['auc_rep_mean'].mean()
sd = auc_by_rep['auc_rep_mean'].std()
median = auc_by_rep['auc_rep_mean'].median()
mn = auc_by_rep['auc_rep_mean'].min(); mx = auc_by_rep['auc_rep_mean'].max()
s2 = [["Statistic","AUC per repeat (n = 100 repeats × 10 outer folds)","Genes selected per outer fold (n = 1000)"],
      ["Mean", f"{mean:.4f}", f"{n_genes['n_genes'].mean():.2f}"],
      ["Standard deviation", f"{sd:.4f}", f"{n_genes['n_genes'].std():.2f}"],
      ["Median", f"{median:.4f}", f"{n_genes['n_genes'].median():.1f}"],
      ["2.5th percentile", f"{q025:.4f}", f"{np.percentile(n_genes['n_genes'],2.5):.0f}"],
      ["97.5th percentile", f"{q975:.4f}", f"{np.percentile(n_genes['n_genes'],97.5):.0f}"],
      ["Minimum", f"{mn:.4f}", f"{int(n_genes['n_genes'].min())}"],
      ["Maximum", f"{mx:.4f}", f"{int(n_genes['n_genes'].max())}"],
      ["Range", f"[{q025:.3f}, {q975:.3f}]",
       f"[{int(n_genes['n_genes'].min())}, {int(n_genes['n_genes'].max())}]"],
      ["Summary line in main text",
       f"AUC = {mean:.3f}, 95 % CI {q025:.3f}–{q975:.3f}",
       f"Median {int(n_genes['n_genes'].median())} genes, range {int(n_genes['n_genes'].min())}–{int(n_genes['n_genes'].max())}"]]
make_three_line_table(doc, s2, font_size=9)

# Table S3
add_caption(doc, "Table S3.  Phenotype summary for each cohort.  Counts are "
                  "based on the analytic samples retained after quality "
                  "control (matches Table 1 totals).")
def sex_summary(df, col='sex'):
    if col not in df.columns: return "n/a"
    vc = df[col].astype(str).str.lower().str[0].value_counts()
    m = int(vc.get('m',0)); f = int(vc.get('f',0)); na = len(df)-m-f
    parts = []
    if m: parts.append(f"M {m}")
    if f: parts.append(f"F {f}")
    if na: parts.append(f"NA {na}")
    return " / ".join(parts) if parts else "n/a"
subtype_notes = {
    'GPL570':'Autistic disorder (clinical) + ASD-PDD-NOS',
    'GPL6244':'Autistic disorder (clinical) + ASD-PDD-NOS',
    'GSE42133':'Toddler ASD (Pierce 2011-2013 cohort)',
    'GSE25507':'ASD diagnosis (Alter cohort)',
}
s3 = [["Cohort","n total","n ASD","n Ctrl","Sex (M/F/NA)","ASD subtype (if reported)"]]
for key, dfp in ph.items():
    n = len(dfp); n_asd = int((dfp['label']==1).sum()); n_ctrl = int((dfp['label']==0).sum())
    s3.append([key, n, n_asd, n_ctrl, sex_summary(dfp), subtype_notes.get(key,'—')])
make_three_line_table(doc, s3, font_size=9)

# Table S4
add_caption(doc, "Table S4.  Full Pearson correlation table of LP₇ (or LP₆ on "
                  "GSE42133) with eight blood-cell-type marker scores in each "
                  "of the four cohorts, including significance (p) and the "
                  "number of markers from the consensus marker list that "
                  "mapped onto the platform.  Marker genes follow the "
                  "consensus CIBERSORT LM22 + xCell list (see Methods).")
cell_types = ['Monocytes','Neutrophils','Dendritic','NK cells',
              'CD4 T cells','CD8 T cells','B cells','Erythrocytes']
canonical_n = {'Monocytes':8,'Neutrophils':8,'Dendritic':4,'NK cells':6,
               'CD4 T cells':7,'CD8 T cells':6,'B cells':8,'Erythrocytes':5}
cohort_keys = [
    ('GSE18123-GPL570', 'GPL570 (whole blood, training)', cm['correlations']),
    ('GSE18123-GPL6244', 'GPL6244 (whole blood, ext)', cm['correlations']),
    ('GSE42133', 'GSE42133', cm42),
    ('GSE25507', 'GSE25507 (lymphocytes)', cm['correlations']),
]
def fmt(r, p):
    rs = f"{r:+.3f}"
    if p < 0.0001: ps = "< 0.0001"
    elif p < 0.001: ps = f"{p:.4f}"
    elif p < 0.01: ps = f"{p:.3f}"
    else: ps = f"{p:.3f}"
    return f"{rs} ; p = {ps}"
s4 = [["Cell type","Markers used (consensus)",
       "GSE18123-GPL570  (r ; p)", "GSE18123-GPL6244  (r ; p)",
       "GSE42133  (r ; p)", "GSE25507  (r ; p)"]]
for ct in cell_types:
    row = [ct, str(canonical_n[ct])]
    for cohort_id, key, src in cohort_keys:
        info = src.get(ct) if cohort_id == 'GSE42133' else src.get(key, {}).get(ct)
        row.append(fmt(info['r'], info['p']) if info else "n/a")
    s4.append(row)
make_three_line_table(doc, s4, font_size=8)

# Table S5
add_caption(doc, "Table S5.  Direction comparison between the locked panel's "
                  "training-cohort coefficients and the within-tissue refit "
                  "coefficients on GSE25507.  Only 2 of 7 genes retain their "
                  "training direction in the lymphocyte-enriched cohort, "
                  "consistent with a partial inversion of the cellular "
                  "substrate the score tracks; this is the basis for our "
                  "framing of GSE25507 as a tissue-boundary exploration "
                  "rather than transportability.")
s5 = [["Gene","Locked panel coefficient (training)",
       "ASD-vs-control Δ in training (log₂)","ASD-vs-control Δ in GSE25507 (log₂)",
       "t (GSE25507)","p (GSE25507)","Direction concordant?"]]
for r in g25['direction_table']:
    s5.append([r['gene'], f"{r['coef_train']:+.4f}",
                f"{r['delta_train']:+.3f}", f"{r['delta_GSE25507']:+.3f}",
                f"{r['t_GSE25507']:+.3f}", f"{r['p_GSE25507']:.3f}",
                "Yes" if r['concordant'] else "No"])
s5.append(["Summary","—","—","—","—","—",
            f"{g25['n_genes_concordant']} / {g25['n_genes_total']} concordant"])
make_three_line_table(doc, s5, font_size=8)

foot = doc.add_paragraph()
foot.paragraph_format.space_before = Pt(8)
fr = foot.add_run("All numerical values were extracted from the locked "
                   "pipeline artifacts shipped with the asd-rerun code "
                   "package, and can be reproduced by running 04_nested_cv.py, "
                   "10_tissue_boundary_GSE25507.py and 11_cell_markers.py "
                   "with seed = 42.")
fr.font.size = Pt(9); fr.font.italic = True; fr.font.name = 'Times New Roman'

doc.save('../Supplementary_Tables.docx')
print("Supplementary tables written")
