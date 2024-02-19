import pandas as pd
import os
import numpy as np

def TCGA(file_path):
    sample_annotations = pd.read_table(file_path, index_col=0, na_values="'--")

    sample_annotations = sample_annotations.loc[:, sample_annotations.isna().mean(axis=0) < .9]
    sample_annotations = sample_annotations.drop(['treatment_or_therapy', 'treatment_type', 
                              'tumor_grade', 'classification_of_tumor', 'last_known_disease_status', 'prior_malignancy', 'prior_treatment', 'progression_or_recurrence', 'site_of_resection_or_biopsy', 'synchronous_malignancy',
                              'case_submitter_id', 'project_id', 'ethnicity', 'gender', 'days_to_birth', 'icd_10_code', 'tissue_or_organ_of_origin'],
                             axis=1)
    sample_annotations = sample_annotations.drop_duplicates()
    assert not sample_annotations.index.duplicated().any()

    sample_annotations['ajcc_pathologic_t'] = sample_annotations['ajcc_pathologic_t'].apply(lambda x:x[:2])
    sample_annotations['ajcc_pathologic_stage'] = sample_annotations['ajcc_pathologic_stage'].str.rstrip('ABC')

    return sample_annotations

