import run_params as rp
from src.data.parse_overview import parse_sample_overview
from src.data.preprocessing import condense_seq_data, add_overview_metrics
from src.data.clustering import cluster_each_patient_seq, \
    cluster_all_patient_seqs, add_cancer_assoc_annot
from src.data.db_annotations import get_deepcat_annotations, annotate_overview
from src.analysis.group_associations import assoc_btwn_groups, \
    assoc_within_groups, patients_in_supercluster
from src.visualization.generate_plots import plot_pop_metrics, \
    plot_cancer_assoc_tcrs
def main():
    ov_path = './data/analysis/SampleOverview.tsv'
    ov_parsed_path = './data/processed/overview_parsed.csv'
    raw_patient_dir = './data/raw/ChemoProjTCRs'
    seq_condensed_dir = './data/processed/condensed_seq_data'
    ov_added_metrics_path = 'data/processed/overview_with_metrics.csv'
    each_pt_cluster_dir = 'data/processed/each_patient_clustered'
    all_pt_cluster_path = 'data/processed/all_patient_clustering.csv'
    deepcat_db_path = 'data/raw/deepcat_cdr3.txt'
    deepcat_annotate_dir = 'data/processed/cancer_annotated'
    ov_num_cancer_tcrs_path = 'data/processed/overview_with_cancer_tcrs.csv'
    overlap_with_supercluster_dir = 'data/processed/overlap_with_supercluster'
    group_assoc_summaries_dir = 'data/processed/group_assoc_summaries'
    btwn_group_summary_dir = 'data/processed/btwn_group_summary'
    
    
    # ---------------- PREPROCESSING ---------------- 
    if rp.preprocess_all or rp.parse_overview:
        print('Parsing overview file...')
        parse_sample_overview(ov_path, ov_parsed_path)
        print(f'    Parsed overview saved to {ov_parsed_path}\n')
    if rp.preprocess_all or rp.condense_seqs:
        print('Condensing sequencing data for each patient...')
        condense_seq_data(raw_patient_dir, seq_condensed_dir)
        print(f'    Condensed patient data saved to {seq_condensed_dir}\n')

    # ---------------- ANALYSIS ---------------- 
    # POPULATION-LEVEL METRICS
    if rp.analyze_all or rp.calc_pop_metrics:
        print('Calculating/saving population metrics for each patient...')
        add_overview_metrics(ov_parsed_path, seq_condensed_dir, 
                             ov_added_metrics_path)
        print('    Overview with population metrics saved to ' + 
              f'{ov_added_metrics_path}\n')
    
    # CLUSTERING
    if rp.analyze_all or rp.cluster_each_patient:
        print("Clustering each patient's data...")
        cluster_each_patient_seq(seq_condensed_dir, each_pt_cluster_dir)
        print('    Clustered patient data (each patient) saved to ' + 
              f'{each_pt_cluster_dir}\n')
    if rp.analyze_all or rp.cluster_all_patients:
        print('Clustering sequences across all patients...')
        # cluster_all_patient_seqs(seq_condensed_dir, all_pt_cluster_path)
        add_cancer_assoc_annot(deepcat_db_path, all_pt_cluster_path)
        print('    Clustered patient data (all patients combined) saved to ' + 
              f'{all_pt_cluster_path}\n')
    
    # DEEPCAT CANCER ASSOCIATION ANNOTATIONS
    if rp.analyze_all or rp.annotate_patients:
        print("Annotating each patient's sequences for cancer association...")
        get_deepcat_annotations(deepcat_db_path, seq_condensed_dir, 
                                deepcat_annotate_dir)
        print('    Clustered patient data (each patient) saved to ' + 
              f'{deepcat_annotate_dir}\n')
        
    # CALCULATE CANCER ASSOCIATED TCRS FOR EACH PATIENT
    if rp.analyze_all or rp.cancer_tcrs_each_pt:
        print('Finding number of cancer-associated TCRs for each patient...')
        annotate_overview(ov_added_metrics_path, deepcat_annotate_dir, 
                         each_pt_cluster_dir, seq_condensed_dir, 
                         ov_num_cancer_tcrs_path)
        print('    Number of cancer-assocated TCRs for each patient saved to ' + 
              f'{ov_num_cancer_tcrs_path}\n')
        
    # FIND ASSOCIATIONS WITHIN/BETWEEN PATIENT GROUPS
    if rp.analyze_all or rp.pts_in_supercluster:
        print('Identifying patient sequences in supercluster...')
        patients_in_supercluster(all_pt_cluster_path, seq_condensed_dir, overlap_with_supercluster_dir)

    if rp.analyze_all or rp.calc_associations_within_groups:
        print('Finding associations within patient groups...')
        assoc_within_groups(ov_num_cancer_tcrs_path, overlap_with_supercluster_dir, group_assoc_summaries_dir)
    if rp.analyze_all or rp.calc_associations_between_groups:
        print('Finding associations between patient groups...')
        assoc_btwn_groups(group_assoc_summaries_dir, btwn_group_summary_dir)


    # ---------------- VISUALIZATION ---------------- 
    if rp.visualize_all or rp.vis_population_metrics:
        print('Plotting population-level metrics...')
        plot_pop_metrics()
    if rp.visualize_all or rp.vis_num_cancer_tcrs:
        print('Plotting number of cancer-associated TCRs per group...')
        plot_cancer_assoc_tcrs()


if __name__ == '__main__':
    main()