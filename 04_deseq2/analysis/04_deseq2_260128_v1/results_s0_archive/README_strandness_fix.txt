変更日: 2026-04-19
変更内容: featureCounts strandness修正 (-s 0 → -s 2)
理由: NEBNext Ultra II Directional (dUTP, RF) ライブラリに対してunstranded設定で
      実行していたことが判明。infer_experiment.py で RF strand-specific (90%)を確認。
旧結果: results_s0_archive/ に保存
