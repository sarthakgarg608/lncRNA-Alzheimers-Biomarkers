# Stage 10 preflight audit

Stage 09 counted both events and individuals. Grouping the saved `outputs/tables/stage09/sample_qc_flags.csv` by dataset and counting rows versus distinct sample_id gives the following totals.

| Cohort | Flag events | Flagged samples |
|---|---:|---:|
| GSE63061 |64|58|
| GSE63060 |20|20|
| GSE140829 |55|47|
| GSE125583 |63|53|
| Total |202|178|

The 63 brain events comprise 41 extreme library sizes and 22 high zero fractions. All sample rows remain eligible for annotation.

The GSE125583 GEO series matrix contains `!Sample_characteristics_ch1` rows keyed by `age` and `Sex`. Stage 08's narrative claim that these fields were absent was incorrect. Script 01 accepts age and age at death as alternative keys and exports sample_id and age. Stage 09 script 05 reads the saved brain metadata.csv, restricts diagnosis to AD/Control, coerces age to numeric, drops missing ages for this test alone and calculates scipy.stats.f_oneway and kruskal. GSM3577605 has the literal age NA. Of 289 samples, 218 AD and 70 controls supplied ages. Means were 83.376 and 87.000 years respectively. Recalculation returned F=11.643883, p=0.000737123 and H=9.761713, p=0.00178184. The published p-value is the ANOVA result. Stage 09 reversed the direction and inferred a Braak association without a corresponding test.

The saved batch_effect_assessment.csv attributes 35.96% of PC1 score variance and 42.72% of PC2 score variance to 19 external batches. Each percentage has its own PC score variance denominator. The batch–diagnosis relationship remains untested. Only the flag and batch tables were needed to audit the external cohort.

