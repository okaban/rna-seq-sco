#!/usr/bin/env python3
"""
H29: Quantitative Boundary Between Shielded and Exposed Regulators
===================================================================

Builds feature matrix for 1,055 regulatory genes, performs ROC analysis,
logistic regression, decision tree classification, expression-distance
scatter with LOWESS, cross-validation, and expression quintile analysis
to determine the boundary between "shielded" and "exposed" regulators.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import to_rgba
import matplotlib.patches as mpatches
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier, export_text
from sklearn.metrics import roc_curve, auc, roc_auc_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from scipy import stats
import statsmodels.api as sm
from statsmodels.nonparametric.smoothers_lowess import lowess

# ============================================================
# Configuration
# ============================================================
BASE = Path('/Users/okaban/bioinfo/rna-seq')
ANALYSIS_DIR = BASE / '11_epigenome_integration/analysis/52_shielded_exposed_boundary'
FIG_DIR = ANALYSIS_DIR / 'figures'
TBL_DIR = ANALYSIS_DIR / 'tables'

CHROM_LEN = 8_667_507
ARM_LEFT_MAX = 1_500_000
ARM_RIGHT_MIN = 7_167_508

SEED = 42
np.random.seed(SEED)

# ============================================================
# 1. Load data
# ============================================================
print("=" * 70)
print("H29: Quantitative Boundary Between Shielded and Exposed Regulators")
print("=" * 70)

# Gene-level metrics from H27
glm = pd.read_csv(
    BASE / '11_epigenome_integration/analysis/50_coordinated_regulators_protection/tables/gene_level_metrics.tsv',
    sep='\t'
)
print(f"\nGene-level metrics: {len(glm)} genes, columns: {list(glm.columns)}")

# 62 coordinated regulators
coord = pd.read_csv(
    BASE / '11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/coordinated_regulatory_genes.tsv',
    sep='\t'
)
coord_tags = set(coord['locus_tag'].tolist())
print(f"Coordinated regulators: {len(coord_tags)} genes")

# All 1,055 regulatory genes
all_reg = pd.read_csv(
    BASE / '11_epigenome_integration/analysis/29_genomewide_TF_screen/tables/all_regulatory_genes.tsv',
    sep='\t'
)
all_reg_tags = set(all_reg['locus_tag'].tolist())
print(f"All regulatory genes: {len(all_reg_tags)} genes")

# DESeq2 results
deseq_t2 = pd.read_csv(
    BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_2_vs_1.tsv',
    sep='\t'
)
deseq_t3 = pd.read_csv(
    BASE / '04_deseq2/analysis/04_deseq2_260128_v1/results/DESeq2_M145_3_vs_1.tsv',
    sep='\t'
)
# Rename for merge
deseq_t2 = deseq_t2.rename(columns={
    'gene_id': 'locus_tag', 'baseMean': 'baseMean',
    'log2FoldChange': 'LFC_T2vsT1', 'padj': 'padj_T2vsT1_deseq2'
})[['locus_tag', 'baseMean', 'LFC_T2vsT1', 'padj_T2vsT1_deseq2']]
deseq_t3 = deseq_t3.rename(columns={
    'gene_id': 'locus_tag',
    'log2FoldChange': 'LFC_T3vsT1', 'padj': 'padj_T3vsT1_deseq2'
})[['locus_tag', 'LFC_T3vsT1', 'padj_T3vsT1_deseq2']]

# Gene annotation for gene_length
annot = pd.read_csv(
    BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv',
    sep='\t'
)
annot = annot.rename(columns={'gene_id': 'locus_tag'})
annot['gene_length'] = annot['end'] - annot['start'] + 1

# FIMO results
fimo_lines = []
with open(BASE / '13_TF_binding-site/analysis/01_master_TF_list_260206_v1/intermediate/fimo_results/fimo.gff') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        start = int(parts[3])
        end = int(parts[4])
        fimo_lines.append({'fimo_start': start, 'fimo_end': end})
fimo_df = pd.DataFrame(fimo_lines)
print(f"FIMO hits loaded: {len(fimo_df)}")

# ============================================================
# 2. Build feature matrix for 1,055 regulatory genes
# ============================================================
print("\n--- Building feature matrix ---")

# Start from gene_level_metrics (which has all 1,055)
df = glm.copy()

# Map exposed vs shielded
df['is_exposed'] = df['locus_tag'].isin(coord_tags).astype(int)
print(f"Exposed: {df['is_exposed'].sum()}, Shielded: {(df['is_exposed'] == 0).sum()}")

# Merge DESeq2 expression data
df = df.merge(deseq_t2, on='locus_tag', how='left')
df = df.merge(deseq_t3, on='locus_tag', how='left')

# Merge gene length from annotation
annot_sub = annot[['locus_tag', 'gene_length']].drop_duplicates(subset='locus_tag')
df = df.merge(annot_sub, on='locus_tag', how='left')

# If gene_length missing, compute from start/end in glm
mask_no_len = df['gene_length'].isna()
if mask_no_len.any():
    df.loc[mask_no_len, 'gene_length'] = df.loc[mask_no_len, 'end'] - df.loc[mask_no_len, 'start'] + 1

# Compute region (core=0, arm=1)
def assign_region(row):
    midpoint = (row['start'] + row['end']) / 2
    if midpoint <= ARM_LEFT_MAX or midpoint >= ARM_RIGHT_MIN:
        return 1
    return 0
df['region'] = df.apply(assign_region, axis=1)

# Log-transform expression
df['log2_baseMean_p1'] = np.log2(df['baseMean'].fillna(0) + 1)

# Rename H27 columns to standard names
df = df.rename(columns={
    'nearest_methyl_dist': 'nearest_methyl_distance',
    'sites_within_2kb': 'n_methyl_sites_2kb'
})

# Count FIMO hits in promoter (-500 to +100 from TSS) for each gene
def count_fimo_promoter(gene_row, fimo_positions):
    tss = gene_row['tss']
    strand = gene_row['strand']
    if strand == '+':
        prom_start = tss - 500
        prom_end = tss + 100
    else:
        prom_start = tss - 100
        prom_end = tss + 500
    # Count FIMO hits overlapping promoter
    hits = fimo_positions[(fimo_positions['fimo_start'] <= prom_end) &
                          (fimo_positions['fimo_end'] >= prom_start)]
    return len(hits)

print("Counting FIMO hits per promoter...")
df['n_FIMO_hits'] = df.apply(lambda r: count_fimo_promoter(r, fimo_df), axis=1)

# Absolute LFC
df['abs_LFC_T2vsT1'] = df['LFC_T2vsT1'].abs()
df['abs_LFC_T3vsT1'] = df['LFC_T3vsT1'].abs()

# Final feature set
feature_cols = [
    'baseMean', 'log2_baseMean_p1', 'nearest_methyl_distance',
    'n_methyl_sites_2kb', 'gene_length', 'region',
    'LFC_T2vsT1', 'LFC_T3vsT1', 'n_FIMO_hits'
]

# Drop rows with essential NAs
essential_cols = ['baseMean', 'nearest_methyl_distance', 'n_methyl_sites_2kb', 'gene_length']
n_before = len(df)
df_analysis = df.dropna(subset=essential_cols).copy()
print(f"After dropping NAs: {len(df_analysis)} / {n_before} genes")
print(f"  Exposed: {df_analysis['is_exposed'].sum()}")
print(f"  Shielded: {(df_analysis['is_exposed'] == 0).sum()}")

# Fill remaining NAs in LFC with 0
df_analysis['LFC_T2vsT1'] = df_analysis['LFC_T2vsT1'].fillna(0)
df_analysis['LFC_T3vsT1'] = df_analysis['LFC_T3vsT1'].fillna(0)

# ============================================================
# 3. Save full feature matrix
# ============================================================
save_cols = ['locus_tag', 'gene_name', 'old_locus_tag', 'product', 'tf_family',
             'start', 'end', 'strand', 'tss', 'region', 'is_exposed',
             'baseMean', 'log2_baseMean_p1', 'nearest_methyl_distance',
             'n_methyl_sites_2kb', 'gene_length', 'LFC_T2vsT1', 'LFC_T3vsT1',
             'n_FIMO_hits']
existing_save_cols = [c for c in save_cols if c in df_analysis.columns]
df_analysis[existing_save_cols].to_csv(TBL_DIR / 'all_genes_features.tsv', sep='\t', index=False)
print(f"\nSaved feature matrix: {TBL_DIR / 'all_genes_features.tsv'}")

# ============================================================
# 4. Univariate ROC Analysis
# ============================================================
print("\n--- Univariate ROC Analysis ---")

y = df_analysis['is_exposed'].values

roc_results = []
roc_curves = {}

# For each feature, calculate ROC (try both directions: higher = exposed, lower = exposed)
for feat in feature_cols:
    vals = df_analysis[feat].values
    # Remove NaN for this feature
    mask = ~np.isnan(vals)
    if mask.sum() < 50:
        continue
    vals_clean = vals[mask]
    y_clean = y[mask]

    # Calculate AUC
    try:
        fpr, tpr, thresholds = roc_curve(y_clean, vals_clean)
        auc_val = auc(fpr, tpr)
    except:
        continue

    # If AUC < 0.5, flip (lower values predict exposed)
    direction = 'higher_exposed'
    if auc_val < 0.5:
        auc_val = 1 - auc_val
        fpr, tpr, thresholds = roc_curve(y_clean, -vals_clean)
        auc_val_check = auc(fpr, tpr)
        direction = 'lower_exposed'

    # Youden index
    youden = tpr - fpr
    best_idx = np.argmax(youden)
    best_threshold = thresholds[best_idx]
    if direction == 'lower_exposed':
        best_threshold = -best_threshold  # flip back
    best_sensitivity = tpr[best_idx]
    best_specificity = 1 - fpr[best_idx]

    # PPV, NPV at optimal threshold
    if direction == 'higher_exposed':
        pred_positive = vals_clean >= best_threshold
    else:
        pred_positive = vals_clean <= best_threshold
    tp = ((pred_positive) & (y_clean == 1)).sum()
    fp = ((pred_positive) & (y_clean == 0)).sum()
    tn = ((~pred_positive) & (y_clean == 0)).sum()
    fn = ((~pred_positive) & (y_clean == 1)).sum()
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0

    # Bootstrap 95% CI for AUC
    n_boot = 1000
    boot_aucs = []
    for _ in range(n_boot):
        idx = np.random.choice(len(y_clean), size=len(y_clean), replace=True)
        y_b = y_clean[idx]
        v_b = vals_clean[idx]
        if len(np.unique(y_b)) < 2:
            continue
        try:
            if direction == 'lower_exposed':
                boot_aucs.append(roc_auc_score(y_b, -v_b))
            else:
                boot_aucs.append(roc_auc_score(y_b, v_b))
        except:
            continue
    ci_lo = np.percentile(boot_aucs, 2.5) if boot_aucs else np.nan
    ci_hi = np.percentile(boot_aucs, 97.5) if boot_aucs else np.nan

    roc_results.append({
        'feature': feat,
        'AUC': round(auc_val, 4),
        'AUC_CI_lo': round(ci_lo, 4),
        'AUC_CI_hi': round(ci_hi, 4),
        'direction': direction,
        'optimal_threshold': round(best_threshold, 4),
        'sensitivity': round(best_sensitivity, 4),
        'specificity': round(best_specificity, 4),
        'PPV': round(ppv, 4),
        'NPV': round(npv, 4),
        'Youden_J': round(best_sensitivity + best_specificity - 1, 4),
        'n_valid': int(mask.sum())
    })

    roc_curves[feat] = (fpr, tpr, auc_val, direction)

roc_df = pd.DataFrame(roc_results).sort_values('AUC', ascending=False)
roc_df.to_csv(TBL_DIR / 'ROC_analysis.tsv', sep='\t', index=False)
print("\nROC Analysis Results:")
print(roc_df.to_string(index=False))

# ============================================================
# 5. Nearest methylation distance distribution
# ============================================================
print("\n--- Nearest methylation distance distribution ---")
exposed = df_analysis[df_analysis['is_exposed'] == 1]['nearest_methyl_distance']
shielded = df_analysis[df_analysis['is_exposed'] == 0]['nearest_methyl_distance']
print(f"Exposed (n={len(exposed)}): median={exposed.median():.0f}, mean={exposed.mean():.0f}")
print(f"Shielded (n={len(shielded)}): median={shielded.median():.0f}, mean={shielded.mean():.0f}")
mwu_stat, mwu_p = stats.mannwhitneyu(exposed, shielded, alternative='two-sided')
print(f"Mann-Whitney U: stat={mwu_stat:.0f}, p={mwu_p:.2e}")

# ============================================================
# 6. Multivariate Logistic Regression
# ============================================================
print("\n--- Multivariate Logistic Regression ---")

lr_features = ['log2_baseMean_p1', 'gene_length', 'region', 'n_FIMO_hits']
X_lr = df_analysis[lr_features].copy()
y_lr = df_analysis['is_exposed'].values

# Fill NAs
X_lr = X_lr.fillna(0)

# Statsmodels for p-values
X_lr_sm = sm.add_constant(X_lr)
try:
    logit_model = sm.Logit(y_lr, X_lr_sm).fit(disp=0, maxiter=100)
    print(logit_model.summary2())

    # Extract results
    lr_results = pd.DataFrame({
        'feature': ['intercept'] + lr_features,
        'coefficient': logit_model.params.values,
        'std_error': logit_model.bse.values,
        'z_value': logit_model.tvalues.values,
        'p_value': logit_model.pvalues.values,
        'odds_ratio': np.exp(logit_model.params.values),
        'OR_CI_lo': np.exp(logit_model.conf_int().iloc[:, 0].values),
        'OR_CI_hi': np.exp(logit_model.conf_int().iloc[:, 1].values)
    })

    # Standardized coefficients for feature importance
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_lr)
    logit_scaled = sm.Logit(y_lr, sm.add_constant(X_scaled)).fit(disp=0, maxiter=100)
    lr_results['standardized_coef'] = [np.nan] + list(logit_scaled.params[1:])

    lr_results.to_csv(TBL_DIR / 'logistic_regression.tsv', sep='\t', index=False)
    print("\nLogistic Regression Results:")
    print(lr_results.to_string(index=False))

    # Multivariate AUC
    y_pred_prob = logit_model.predict(X_lr_sm)
    mv_auc = roc_auc_score(y_lr, y_pred_prob)
    print(f"\nMultivariate AUC: {mv_auc:.4f}")

    # Store for plotting
    fpr_mv, tpr_mv, _ = roc_curve(y_lr, y_pred_prob)

except Exception as e:
    print(f"Logistic regression error: {e}")
    lr_results = None
    mv_auc = None
    fpr_mv, tpr_mv = None, None

# ============================================================
# 7. Decision Tree
# ============================================================
print("\n--- Decision Tree (max_depth=3) ---")

dt_features = ['baseMean', 'log2_baseMean_p1', 'nearest_methyl_distance',
               'n_methyl_sites_2kb', 'gene_length', 'region', 'n_FIMO_hits']
X_dt = df_analysis[dt_features].fillna(0).values
y_dt = y.copy()

# Fit trees at depth 2 and 3
for depth in [2, 3]:
    dt = DecisionTreeClassifier(max_depth=depth, random_state=SEED, class_weight='balanced')
    dt.fit(X_dt, y_dt)
    dt_auc = roc_auc_score(y_dt, dt.predict_proba(X_dt)[:, 1])
    print(f"\nDecision Tree (depth={depth}), Training AUC: {dt_auc:.4f}")
    tree_text = export_text(dt, feature_names=dt_features, max_depth=depth)
    print(tree_text)

# Use depth=3 for final
dt_final = DecisionTreeClassifier(max_depth=3, random_state=SEED, class_weight='balanced')
dt_final.fit(X_dt, y_dt)
dt_proba = dt_final.predict_proba(X_dt)[:, 1]
dt_auc_final = roc_auc_score(y_dt, dt_proba)

# Feature importances
dt_importances = pd.DataFrame({
    'feature': dt_features,
    'importance': dt_final.feature_importances_
}).sort_values('importance', ascending=False)
print("\nDecision Tree Feature Importances:")
print(dt_importances.to_string(index=False))

# ============================================================
# 8. Expression vs Distance Scatter with LOWESS
# ============================================================
print("\n--- Expression vs Distance scatter ---")

x_expr = df_analysis['baseMean'].values
y_dist = df_analysis['nearest_methyl_distance'].values
labels = df_analysis['is_exposed'].values

# LOWESS fit (sort by x for smooth curve)
sort_idx = np.argsort(x_expr)
x_sorted = x_expr[sort_idx]
y_sorted = y_dist[sort_idx]

# Use log scale for better LOWESS
x_log = np.log10(x_sorted + 1)
lowess_result = lowess(y_sorted, x_log, frac=0.3, return_sorted=True)

# Spearman correlation
rho, p_spearman = stats.spearmanr(x_expr, y_dist)
print(f"Spearman rho (baseMean vs nearest_methyl_distance): {rho:.4f}, p={p_spearman:.2e}")

rho_log, p_log = stats.spearmanr(np.log2(x_expr + 1), y_dist)
print(f"Spearman rho (log2(baseMean+1) vs nearest_methyl_distance): {rho_log:.4f}, p={p_log:.2e}")

# ============================================================
# 9. Cross-validation
# ============================================================
print("\n--- 5-fold Cross-Validation ---")

cv_results = []
models_cv = {
    'Logistic (4 feat)': ('logistic', lr_features),
    'Logistic (all)': ('logistic', [f for f in feature_cols if f not in ['baseMean']]),  # exclude baseMean since we have log version
    'Decision Tree (d=3)': ('tree', dt_features),
    'Single: nearest_methyl_distance': ('single', ['nearest_methyl_distance']),
    'Single: n_methyl_sites_2kb': ('single', ['n_methyl_sites_2kb']),
    'Single: baseMean': ('single', ['baseMean']),
}

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)

for model_name, (model_type, feats) in models_cv.items():
    X_cv = df_analysis[feats].fillna(0).values
    y_cv = y.copy()

    fold_aucs = []
    for train_idx, test_idx in skf.split(X_cv, y_cv):
        X_train, X_test = X_cv[train_idx], X_cv[test_idx]
        y_train, y_test = y_cv[train_idx], y_cv[test_idx]

        if model_type == 'logistic':
            clf = LogisticRegression(max_iter=1000, random_state=SEED, class_weight='balanced')
            clf.fit(X_train, y_train)
            prob = clf.predict_proba(X_test)[:, 1]
        elif model_type == 'tree':
            clf = DecisionTreeClassifier(max_depth=3, random_state=SEED, class_weight='balanced')
            clf.fit(X_train, y_train)
            prob = clf.predict_proba(X_test)[:, 1]
        elif model_type == 'single':
            # For single feature, use direction from ROC analysis
            feat_name = feats[0]
            roc_row = roc_df[roc_df['feature'] == feat_name]
            if len(roc_row) > 0 and roc_row.iloc[0]['direction'] == 'lower_exposed':
                prob = -X_test.ravel()
            else:
                prob = X_test.ravel()

        try:
            fold_auc = roc_auc_score(y_test, prob)
            fold_aucs.append(fold_auc)
        except:
            pass

    mean_auc = np.mean(fold_aucs) if fold_aucs else np.nan
    std_auc = np.std(fold_aucs) if fold_aucs else np.nan

    cv_results.append({
        'model': model_name,
        'n_features': len(feats),
        'mean_AUC': round(mean_auc, 4),
        'std_AUC': round(std_auc, 4),
        'min_AUC': round(min(fold_aucs), 4) if fold_aucs else np.nan,
        'max_AUC': round(max(fold_aucs), 4) if fold_aucs else np.nan,
    })

cv_df = pd.DataFrame(cv_results).sort_values('mean_AUC', ascending=False)
cv_df.to_csv(TBL_DIR / 'cross_validation.tsv', sep='\t', index=False)
print("\nCross-Validation Results:")
print(cv_df.to_string(index=False))

# ============================================================
# 10. Expression Quintile Analysis
# ============================================================
print("\n--- Expression Quintile Analysis ---")

df_analysis['baseMean_quintile'] = pd.qcut(df_analysis['baseMean'], q=5, labels=False, duplicates='drop') + 1

quintile_results = []
for q in sorted(df_analysis['baseMean_quintile'].unique()):
    sub = df_analysis[df_analysis['baseMean_quintile'] == q]
    n_total = len(sub)
    n_exposed = sub['is_exposed'].sum()
    frac_exposed = n_exposed / n_total if n_total > 0 else 0
    mean_dist = sub['nearest_methyl_distance'].mean()
    median_dist = sub['nearest_methyl_distance'].median()
    mean_expr = sub['baseMean'].mean()
    median_expr = sub['baseMean'].median()
    mean_sites = sub['n_methyl_sites_2kb'].mean()

    quintile_results.append({
        'quintile': q,
        'n_genes': n_total,
        'n_exposed': n_exposed,
        'frac_exposed': round(frac_exposed, 4),
        'mean_baseMean': round(mean_expr, 1),
        'median_baseMean': round(median_expr, 1),
        'mean_nearest_methyl_distance': round(mean_dist, 1),
        'median_nearest_methyl_distance': round(median_dist, 1),
        'mean_n_methyl_sites_2kb': round(mean_sites, 2)
    })

quint_df = pd.DataFrame(quintile_results)
quint_df.to_csv(TBL_DIR / 'expression_quintile.tsv', sep='\t', index=False)
print("\nExpression Quintile Analysis:")
print(quint_df.to_string(index=False))

# Jonckheere-Terpstra trend test (approximation with Spearman on quintile vs frac_exposed)
# More rigorous: use the ordered quintile assignments to test trend
from itertools import combinations
def jonckheere_terpstra(groups, values):
    """Manual JT test: count concordant pairs across ordered groups."""
    unique_groups = sorted(set(groups))
    jt_stat = 0
    n_pairs = 0
    for i, j in combinations(range(len(unique_groups)), 2):
        g_i = unique_groups[i]
        g_j = unique_groups[j]
        vals_i = values[groups == g_i]
        vals_j = values[groups == g_j]
        for vi in vals_i:
            for vj in vals_j:
                if vj > vi:
                    jt_stat += 1
                elif vj == vi:
                    jt_stat += 0.5
                n_pairs += 1
    # Z approximation
    n = len(values)
    group_sizes = [np.sum(groups == g) for g in unique_groups]
    E_JT = (n**2 - sum(ni**2 for ni in group_sizes)) / 4
    var_num = n**2 * (2*n + 3)
    for ni in group_sizes:
        var_num -= ni**2 * (2*ni + 3)
    var_JT = var_num / 72
    z = (jt_stat - E_JT) / np.sqrt(var_JT) if var_JT > 0 else 0
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return jt_stat, z, p

# Test: exposed fraction trend across quintiles
groups_q = df_analysis['baseMean_quintile'].values
exposed_vals = df_analysis['is_exposed'].values

jt_stat, jt_z, jt_p = jonckheere_terpstra(groups_q, exposed_vals)
print(f"\nJonckheere-Terpstra test (exposed fraction across expression quintiles):")
print(f"  JT stat={jt_stat:.1f}, z={jt_z:.3f}, p={jt_p:.4e}")

# Also test distance trend
jt_stat_d, jt_z_d, jt_p_d = jonckheere_terpstra(groups_q, df_analysis['nearest_methyl_distance'].values)
print(f"Jonckheere-Terpstra test (nearest distance across expression quintiles):")
print(f"  JT stat={jt_stat_d:.1f}, z={jt_z_d:.3f}, p={jt_p_d:.4e}")

# ============================================================
# VISUALIZATIONS
# ============================================================
print("\n--- Generating Figures ---")

# Color scheme
COL_EXPOSED = '#E74C3C'   # red
COL_SHIELDED = '#3498DB'  # blue
COL_LOESS = '#2C3E50'     # dark

# ---- Panel A: ROC Curves ----
fig_roc, ax_roc = plt.subplots(figsize=(8, 7))

# Top features by AUC
top_features = roc_df.head(6)['feature'].tolist()
colors_roc = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12', '#9B59B6', '#1ABC9C']

for i, feat in enumerate(top_features):
    if feat in roc_curves:
        fpr_f, tpr_f, auc_f, dir_f = roc_curves[feat]
        label = f"{feat} (AUC={auc_f:.3f})"
        ax_roc.plot(fpr_f, tpr_f, color=colors_roc[i % len(colors_roc)],
                    lw=2, label=label)

# Add multivariate
if fpr_mv is not None:
    ax_roc.plot(fpr_mv, tpr_mv, color='black', lw=2.5, ls='--',
                label=f'Multivariate LR (AUC={mv_auc:.3f})')

ax_roc.plot([0, 1], [0, 1], 'k--', alpha=0.3, lw=1)
ax_roc.set_xlabel('False Positive Rate', fontsize=13)
ax_roc.set_ylabel('True Positive Rate (Sensitivity)', fontsize=13)
ax_roc.set_title('ROC Curves: Predicting Exposed vs Shielded Regulators', fontsize=14)
ax_roc.legend(loc='lower right', fontsize=9)
ax_roc.set_xlim([-0.02, 1.02])
ax_roc.set_ylim([-0.02, 1.02])
ax_roc.set_aspect('equal')
fig_roc.tight_layout()
fig_roc.savefig(FIG_DIR / 'ROC_curves.pdf', dpi=300, bbox_inches='tight')
fig_roc.savefig(FIG_DIR / 'ROC_curves.svg', bbox_inches='tight')
plt.close(fig_roc)
print("  Saved ROC_curves.pdf/svg")

# ---- Panel B: Expression vs Distance Scatter ----
fig_scatter, ax_sc = plt.subplots(figsize=(10, 7))

mask_sh = df_analysis['is_exposed'] == 0
mask_ex = df_analysis['is_exposed'] == 1

ax_sc.scatter(df_analysis.loc[mask_sh, 'baseMean'],
              df_analysis.loc[mask_sh, 'nearest_methyl_distance'],
              c=COL_SHIELDED, alpha=0.35, s=20, label=f'Shielded (n={mask_sh.sum()})',
              edgecolors='none', zorder=2)
ax_sc.scatter(df_analysis.loc[mask_ex, 'baseMean'],
              df_analysis.loc[mask_ex, 'nearest_methyl_distance'],
              c=COL_EXPOSED, alpha=0.8, s=50, marker='D',
              label=f'Exposed (n={mask_ex.sum()})', edgecolors='black',
              linewidth=0.5, zorder=3)

# LOWESS
ax_sc.plot(10**lowess_result[:, 0] - 1, lowess_result[:, 1],
           color=COL_LOESS, lw=3, label='LOWESS smoothing', zorder=4)

# Optimal threshold line for nearest_methyl_distance
dist_row = roc_df[roc_df['feature'] == 'nearest_methyl_distance']
if len(dist_row) > 0:
    thresh_dist = dist_row.iloc[0]['optimal_threshold']
    ax_sc.axhline(y=thresh_dist, color='gray', ls='--', lw=1.5, alpha=0.7,
                  label=f'Distance threshold: {thresh_dist:.0f} bp')

ax_sc.set_xscale('log')
ax_sc.set_xlabel('baseMean (average expression)', fontsize=13)
ax_sc.set_ylabel('Nearest methylation distance (bp)', fontsize=13)
ax_sc.set_title('Expression vs Methylation Distance\n(1,055 Regulatory Genes)', fontsize=14)
ax_sc.legend(fontsize=10, loc='upper right')

# Add Spearman annotation
ax_sc.text(0.02, 0.02, f'Spearman rho = {rho:.3f} (p = {p_spearman:.2e})',
           transform=ax_sc.transAxes, fontsize=10,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

fig_scatter.tight_layout()
fig_scatter.savefig(FIG_DIR / 'expression_vs_distance_scatter.pdf', dpi=300, bbox_inches='tight')
fig_scatter.savefig(FIG_DIR / 'expression_vs_distance_scatter.svg', bbox_inches='tight')
plt.close(fig_scatter)
print("  Saved expression_vs_distance_scatter.pdf/svg")

# ---- Panel C: Expression Quintile Analysis ----
fig_quint, (ax_q1, ax_q2) = plt.subplots(1, 2, figsize=(12, 5))

x_q = quint_df['quintile'].values

# Left: fraction exposed
bars = ax_q1.bar(x_q, quint_df['frac_exposed'].values * 100, color=COL_EXPOSED, alpha=0.7,
                 edgecolor='darkred', linewidth=1)
for i, (xv, yv) in enumerate(zip(x_q, quint_df['frac_exposed'].values * 100)):
    ax_q1.text(xv, yv + 0.3, f"{yv:.1f}%\n(n={quint_df.iloc[i]['n_exposed']})",
               ha='center', fontsize=9)
ax_q1.set_xlabel('Expression Quintile', fontsize=12)
ax_q1.set_ylabel('% Exposed Regulators', fontsize=12)
ax_q1.set_title('Fraction Exposed by Expression Quintile', fontsize=13)
ax_q1.set_xticks(x_q)
ax_q1.set_xticklabels([f'Q{q}\n({quint_df.iloc[q-1]["mean_baseMean"]:.0f})' for q in x_q], fontsize=9)
ax_q1.text(0.02, 0.95, f'JT z={jt_z:.2f}, p={jt_p:.3e}',
           transform=ax_q1.transAxes, fontsize=9, va='top',
           bbox=dict(boxstyle='round', facecolor='lightyellow'))

# Right: mean nearest distance
ax_q2.bar(x_q, quint_df['mean_nearest_methyl_distance'].values, color=COL_SHIELDED,
          alpha=0.7, edgecolor='navy', linewidth=1)
for i, (xv, yv) in enumerate(zip(x_q, quint_df['mean_nearest_methyl_distance'].values)):
    ax_q2.text(xv, yv + 10, f"{yv:.0f}", ha='center', fontsize=9)
ax_q2.set_xlabel('Expression Quintile', fontsize=12)
ax_q2.set_ylabel('Mean Nearest Methylation Distance (bp)', fontsize=12)
ax_q2.set_title('Methylation Distance by Expression Quintile', fontsize=13)
ax_q2.set_xticks(x_q)
ax_q2.set_xticklabels([f'Q{q}' for q in x_q], fontsize=10)
ax_q2.text(0.02, 0.95, f'JT z={jt_z_d:.2f}, p={jt_p_d:.3e}',
           transform=ax_q2.transAxes, fontsize=9, va='top',
           bbox=dict(boxstyle='round', facecolor='lightyellow'))

fig_quint.tight_layout()
fig_quint.savefig(FIG_DIR / 'expression_quintile_analysis.pdf', dpi=300, bbox_inches='tight')
fig_quint.savefig(FIG_DIR / 'expression_quintile_analysis.svg', bbox_inches='tight')
plt.close(fig_quint)
print("  Saved expression_quintile_analysis.pdf/svg")

# ---- Panel E: Feature Importance ----
fig_feat, ax_feat = plt.subplots(figsize=(8, 5))

if lr_results is not None:
    # Plot odds ratios from logistic regression (excluding intercept)
    lr_plot = lr_results[lr_results['feature'] != 'intercept'].copy()
    lr_plot = lr_plot.sort_values('odds_ratio', ascending=True)

    y_pos = range(len(lr_plot))
    ax_feat.barh(y_pos, lr_plot['odds_ratio'].values, color=COL_SHIELDED, alpha=0.7,
                 edgecolor='navy')
    ax_feat.set_yticks(y_pos)
    ax_feat.set_yticklabels(lr_plot['feature'].values, fontsize=11)
    ax_feat.axvline(x=1.0, color='red', ls='--', lw=1.5, alpha=0.7)
    ax_feat.set_xlabel('Odds Ratio (Exposed vs Shielded)', fontsize=12)
    ax_feat.set_title('Logistic Regression: Feature Odds Ratios', fontsize=13)

    # Add significance markers
    for i, (_, row) in enumerate(lr_plot.iterrows()):
        sig = ''
        if row['p_value'] < 0.001:
            sig = '***'
        elif row['p_value'] < 0.01:
            sig = '**'
        elif row['p_value'] < 0.05:
            sig = '*'
        ax_feat.text(row['odds_ratio'] + 0.02, i,
                     f"OR={row['odds_ratio']:.2f} {sig}\n(p={row['p_value']:.3e})",
                     va='center', fontsize=8)

    # Error bars
    ax_feat.errorbar(lr_plot['odds_ratio'].values, y_pos,
                     xerr=[lr_plot['odds_ratio'].values - lr_plot['OR_CI_lo'].values,
                           lr_plot['OR_CI_hi'].values - lr_plot['odds_ratio'].values],
                     fmt='none', color='black', capsize=3)

fig_feat.tight_layout()
fig_feat.savefig(FIG_DIR / 'feature_importance.pdf', dpi=300, bbox_inches='tight')
fig_feat.savefig(FIG_DIR / 'feature_importance.svg', bbox_inches='tight')
plt.close(fig_feat)
print("  Saved feature_importance.pdf/svg")

# ---- Panel D: Decision Tree Visualization (simplified rule diagram) ----
fig_dt, ax_dt = plt.subplots(figsize=(12, 7))
from sklearn.tree import plot_tree
plot_tree(dt_final, feature_names=dt_features, class_names=['Shielded', 'Exposed'],
          filled=True, rounded=True, ax=ax_dt, fontsize=8, proportion=True,
          impurity=True)
ax_dt.set_title(f'Decision Tree (depth=3, AUC={dt_auc_final:.3f})', fontsize=14)
fig_dt.tight_layout()
fig_dt.savefig(FIG_DIR / 'decision_tree.pdf', dpi=300, bbox_inches='tight')
fig_dt.savefig(FIG_DIR / 'decision_tree.svg', bbox_inches='tight')
plt.close(fig_dt)
print("  Saved decision_tree.pdf/svg")

# ---- Comprehensive Summary Figure ----
fig_comp = plt.figure(figsize=(20, 16))
gs = GridSpec(3, 3, figure=fig_comp, hspace=0.35, wspace=0.35)

# A: ROC curves
ax_a = fig_comp.add_subplot(gs[0, 0:2])
for i, feat in enumerate(top_features[:4]):
    if feat in roc_curves:
        fpr_f, tpr_f, auc_f, dir_f = roc_curves[feat]
        ax_a.plot(fpr_f, tpr_f, color=colors_roc[i], lw=2,
                  label=f"{feat} ({auc_f:.3f})")
if fpr_mv is not None:
    ax_a.plot(fpr_mv, tpr_mv, 'k--', lw=2.5, label=f'Multivar LR ({mv_auc:.3f})')
ax_a.plot([0, 1], [0, 1], 'k--', alpha=0.3, lw=1)
ax_a.set_xlabel('FPR')
ax_a.set_ylabel('TPR')
ax_a.set_title('A. ROC Curves', fontsize=14, fontweight='bold')
ax_a.legend(fontsize=8, loc='lower right')
ax_a.set_xlim([-0.02, 1.02])
ax_a.set_ylim([-0.02, 1.02])
ax_a.set_aspect('equal')

# B: Expression vs Distance
ax_b = fig_comp.add_subplot(gs[0, 2])
ax_b.scatter(df_analysis.loc[mask_sh, 'baseMean'],
             df_analysis.loc[mask_sh, 'nearest_methyl_distance'],
             c=COL_SHIELDED, alpha=0.3, s=15, edgecolors='none')
ax_b.scatter(df_analysis.loc[mask_ex, 'baseMean'],
             df_analysis.loc[mask_ex, 'nearest_methyl_distance'],
             c=COL_EXPOSED, alpha=0.8, s=35, marker='D', edgecolors='black', lw=0.3)
ax_b.plot(10**lowess_result[:, 0] - 1, lowess_result[:, 1],
          color=COL_LOESS, lw=2.5)
ax_b.set_xscale('log')
ax_b.set_xlabel('baseMean')
ax_b.set_ylabel('Nearest methyl dist (bp)')
ax_b.set_title('B. Expression vs Distance', fontsize=14, fontweight='bold')
ax_b.text(0.02, 0.02, f'rho={rho:.3f}', transform=ax_b.transAxes, fontsize=9,
          bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

# C: Quintile - fraction exposed
ax_c = fig_comp.add_subplot(gs[1, 0])
ax_c.bar(x_q, quint_df['frac_exposed'].values * 100, color=COL_EXPOSED, alpha=0.7,
         edgecolor='darkred')
for i, (xv, yv) in enumerate(zip(x_q, quint_df['frac_exposed'].values * 100)):
    ax_c.text(xv, yv + 0.2, f"{yv:.1f}%", ha='center', fontsize=8)
ax_c.set_xlabel('Expression Quintile')
ax_c.set_ylabel('% Exposed')
ax_c.set_title('C. Exposed Fraction by Quintile', fontsize=14, fontweight='bold')
ax_c.set_xticks(x_q)
ax_c.text(0.02, 0.95, f'JT p={jt_p:.2e}', transform=ax_c.transAxes, fontsize=9, va='top',
          bbox=dict(boxstyle='round', facecolor='lightyellow'))

# D: Quintile - distance
ax_d = fig_comp.add_subplot(gs[1, 1])
ax_d.bar(x_q, quint_df['mean_nearest_methyl_distance'].values, color=COL_SHIELDED,
         alpha=0.7, edgecolor='navy')
for i, (xv, yv) in enumerate(zip(x_q, quint_df['mean_nearest_methyl_distance'].values)):
    ax_d.text(xv, yv + 5, f"{yv:.0f}", ha='center', fontsize=8)
ax_d.set_xlabel('Expression Quintile')
ax_d.set_ylabel('Mean Nearest Distance (bp)')
ax_d.set_title('D. Distance by Quintile', fontsize=14, fontweight='bold')
ax_d.set_xticks(x_q)
ax_d.text(0.02, 0.95, f'JT p={jt_p_d:.2e}', transform=ax_d.transAxes, fontsize=9, va='top',
          bbox=dict(boxstyle='round', facecolor='lightyellow'))

# E: Feature importance (OR)
ax_e = fig_comp.add_subplot(gs[1, 2])
if lr_results is not None:
    lr_plot = lr_results[lr_results['feature'] != 'intercept'].sort_values('odds_ratio', ascending=True)
    y_pos = range(len(lr_plot))
    ax_e.barh(y_pos, lr_plot['odds_ratio'].values, color=COL_SHIELDED, alpha=0.7)
    ax_e.set_yticks(y_pos)
    ax_e.set_yticklabels(lr_plot['feature'].values, fontsize=9)
    ax_e.axvline(x=1.0, color='red', ls='--', lw=1.5, alpha=0.7)
    ax_e.errorbar(lr_plot['odds_ratio'].values, y_pos,
                  xerr=[lr_plot['odds_ratio'].values - lr_plot['OR_CI_lo'].values,
                        lr_plot['OR_CI_hi'].values - lr_plot['odds_ratio'].values],
                  fmt='none', color='black', capsize=3)
    for i, (_, row) in enumerate(lr_plot.iterrows()):
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
        ax_e.text(max(row['odds_ratio'], row['OR_CI_hi']) + 0.02, i,
                  f"{sig}", va='center', fontsize=8)
ax_e.set_xlabel('Odds Ratio')
ax_e.set_title('E. Feature Importance (LR)', fontsize=14, fontweight='bold')

# F: Distribution comparison (violin/histogram)
ax_f = fig_comp.add_subplot(gs[2, 0:2])

# Side-by-side histograms of nearest methylation distance
bins = np.linspace(0, 5000, 50)
ax_f.hist(shielded.clip(upper=5000), bins=bins, color=COL_SHIELDED, alpha=0.5,
          label=f'Shielded (n={len(shielded)}, med={shielded.median():.0f})',
          density=True, edgecolor='white')
ax_f.hist(exposed.clip(upper=5000), bins=bins, color=COL_EXPOSED, alpha=0.6,
          label=f'Exposed (n={len(exposed)}, med={exposed.median():.0f})',
          density=True, edgecolor='white')
if len(dist_row) > 0:
    ax_f.axvline(x=thresh_dist, color='black', ls='--', lw=2,
                 label=f'Optimal threshold: {thresh_dist:.0f} bp')
ax_f.set_xlabel('Nearest Methylation Distance (bp)')
ax_f.set_ylabel('Density')
ax_f.set_title('F. Distance Distribution: Exposed vs Shielded', fontsize=14, fontweight='bold')
ax_f.legend(fontsize=9)
ax_f.text(0.02, 0.95, f'Mann-Whitney p={mwu_p:.2e}', transform=ax_f.transAxes,
          fontsize=9, va='top', bbox=dict(boxstyle='round', facecolor='lightyellow'))

# G: Cross-validation summary
ax_g = fig_comp.add_subplot(gs[2, 2])
cv_sorted = cv_df.sort_values('mean_AUC', ascending=True)
y_cv_pos = range(len(cv_sorted))
colors_cv = ['#E74C3C' if 'Single' in m else '#3498DB' if 'Logistic' in m else '#2ECC71'
             for m in cv_sorted['model'].values]
ax_g.barh(y_cv_pos, cv_sorted['mean_AUC'].values, color=colors_cv, alpha=0.7, edgecolor='gray')
ax_g.errorbar(cv_sorted['mean_AUC'].values, y_cv_pos,
              xerr=cv_sorted['std_AUC'].values, fmt='none', color='black', capsize=3)
ax_g.set_yticks(y_cv_pos)
ax_g.set_yticklabels(cv_sorted['model'].values, fontsize=8)
ax_g.axvline(x=0.5, color='gray', ls=':', lw=1)
ax_g.axvline(x=0.7, color='orange', ls='--', lw=1, alpha=0.7, label='AUC=0.7')
ax_g.axvline(x=0.8, color='red', ls='--', lw=1, alpha=0.7, label='AUC=0.8')
for i, (_, row) in enumerate(cv_sorted.iterrows()):
    ax_g.text(row['mean_AUC'] + 0.01, i, f"{row['mean_AUC']:.3f}", va='center', fontsize=8)
ax_g.set_xlabel('5-fold CV AUC')
ax_g.set_title('G. Cross-Validation AUC', fontsize=14, fontweight='bold')
ax_g.legend(fontsize=8, loc='lower right')
ax_g.set_xlim([0.4, 1.0])

fig_comp.suptitle('H29: Quantitative Boundary Between Shielded and Exposed Regulators',
                  fontsize=16, fontweight='bold', y=0.98)
fig_comp.savefig(FIG_DIR / 'H29_comprehensive_summary.pdf', dpi=300, bbox_inches='tight')
fig_comp.savefig(FIG_DIR / 'H29_comprehensive_summary.svg', bbox_inches='tight')
plt.close(fig_comp)
print("  Saved H29_comprehensive_summary.pdf/svg")

# ============================================================
# Final Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

# Best univariate feature
best_feat = roc_df.iloc[0]
print(f"\nBest univariate predictor: {best_feat['feature']}")
print(f"  AUC = {best_feat['AUC']:.4f} (95% CI: {best_feat['AUC_CI_lo']:.4f}-{best_feat['AUC_CI_hi']:.4f})")
print(f"  Optimal threshold: {best_feat['optimal_threshold']:.2f}")
print(f"  Sensitivity: {best_feat['sensitivity']:.3f}, Specificity: {best_feat['specificity']:.3f}")
print(f"  PPV: {best_feat['PPV']:.3f}, NPV: {best_feat['NPV']:.3f}")
print(f"  Youden's J: {best_feat['Youden_J']:.3f}")
print(f"  Direction: {best_feat['direction']}")

if mv_auc is not None:
    print(f"\nMultivariate LR AUC: {mv_auc:.4f}")
    print(f"Improvement over best univariate: {mv_auc - best_feat['AUC']:.4f}")

print(f"\nDecision Tree AUC (depth=3): {dt_auc_final:.4f}")

# Best CV model
best_cv = cv_df.iloc[0]
print(f"\nBest cross-validated model: {best_cv['model']}")
print(f"  Mean CV AUC = {best_cv['mean_AUC']:.4f} +/- {best_cv['std_AUC']:.4f}")

# Expression vs distance relationship
print(f"\nExpression-Distance relationship:")
print(f"  Spearman rho = {rho:.4f} (p = {p_spearman:.2e})")

# Quintile trend
print(f"\nExpression Quintile trend:")
print(f"  Exposed fraction JT: z={jt_z:.3f}, p={jt_p:.4e}")
print(f"  Distance JT: z={jt_z_d:.3f}, p={jt_p_d:.4e}")

# Interpretation
print("\n--- Interpretation ---")
top_auc = best_feat['AUC']
if top_auc > 0.8:
    print("AUC > 0.8: STRONG predictive power - clear boundary exists")
elif top_auc > 0.7:
    print("AUC 0.7-0.8: MODERATE predictive power")
elif top_auc > 0.6:
    print("AUC 0.6-0.7: WEAK predictive power, multiple factors contribute")
else:
    print("AUC < 0.6: NO predictive power for this feature")

if mv_auc is not None and mv_auc > top_auc + 0.05:
    print("Multivariate >> Univariate: Multiple factors needed for prediction")
elif mv_auc is not None:
    print("Multivariate ~ Univariate: Single feature captures most information")

print("\n" + "=" * 70)
print("Output files:")
print(f"  {TBL_DIR / 'ROC_analysis.tsv'}")
print(f"  {TBL_DIR / 'logistic_regression.tsv'}")
print(f"  {TBL_DIR / 'expression_quintile.tsv'}")
print(f"  {TBL_DIR / 'all_genes_features.tsv'}")
print(f"  {TBL_DIR / 'cross_validation.tsv'}")
print(f"  {FIG_DIR / 'ROC_curves.pdf/svg'}")
print(f"  {FIG_DIR / 'expression_vs_distance_scatter.pdf/svg'}")
print(f"  {FIG_DIR / 'expression_quintile_analysis.pdf/svg'}")
print(f"  {FIG_DIR / 'feature_importance.pdf/svg'}")
print(f"  {FIG_DIR / 'decision_tree.pdf/svg'}")
print(f"  {FIG_DIR / 'H29_comprehensive_summary.pdf/svg'}")
print("=" * 70)
