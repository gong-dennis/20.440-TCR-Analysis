from os.path import join, splitext
from os import listdir
import pandas as pd
from tqdm import tqdm


def patients_in_supercluster(all_pt_cluster_path, seq_condensed_dir, 
                             overlap_with_supercluster_dir):
    all_pt_clusters_df = pd.read_csv(all_pt_cluster_path, sep=',', 
                                     low_memory=False)
    aa_to_cluster_dict = dict(zip(all_pt_clusters_df.junction_aa, 
                                  all_pt_clusters_df.cluster))
    aa_to_cancer_assoc_dict = dict(zip(all_pt_clusters_df.junction_aa, 
                                   all_pt_clusters_df.is_cancer_assoc))
    
    pt_seq_names = listdir(seq_condensed_dir)
    for f in tqdm(pt_seq_names):
        pt_info_dict = {}
        pt_seq_file_path = join(seq_condensed_dir, f)
        seq_df = pd.read_csv(pt_seq_file_path, sep='\t', low_memory=False)
        for i, row in seq_df.iterrows():
            aa_seq = row['aminoAcid']
            mult = row['count (templates/reads)']
            if aa_seq in aa_to_cluster_dict.keys():
                clust_num = aa_to_cluster_dict[aa_seq]
                is_cancer_assoc = aa_to_cancer_assoc_dict[aa_seq]
                if clust_num in pt_info_dict.keys():
                    pt_info_dict[clust_num][0] += 1
                    pt_info_dict[clust_num][1] += mult
                else:
                    pt_info_dict[clust_num] = [1, mult, is_cancer_assoc]
        
        df_dict = {'cluster_number': [], 'number_tcrs_in_cluster': [], 
                   'total_multiplicity': [], 'is_cancer_assoc': []}
        for c_num in pt_info_dict:
            df_dict['cluster_number'].append(c_num)
            df_dict['number_tcrs_in_cluster'].append(pt_info_dict[c_num][0])
            df_dict['total_multiplicity'].append(pt_info_dict[c_num][1])
            df_dict['is_cancer_assoc'].append(pt_info_dict[c_num][2])
        data_df = pd.DataFrame(df_dict)
        data_df.to_csv(join(overlap_with_supercluster_dir, f), index=False, 
                        sep='\t')


def assoc_within_groups(ov_num_cancer_tcrs_path, overlap_with_supercluster_dir, 
                        group_assoc_summaries_dir):
    ov_df = pd.read_csv(ov_num_cancer_tcrs_path, sep='\t', low_memory=False)
    all_group_labels = list(pd.unique(ov_df['group_label']))
    labels_reformat = [s.lower().replace(' ', '_') for s in all_group_labels]
    for i, lb in enumerate(all_group_labels):
        temp_dict = {}
        lab_slice = ov_df.loc[ov_df['group_label'] == lb]
        for name in lab_slice['sample_name']:
            path_to_df = join(overlap_with_supercluster_dir, name + '.csv')
            sc_overlap_df = pd.read_csv(path_to_df, sep='\t', low_memory=False)
            for _, row in sc_overlap_df.iterrows():
                clust_num = row['cluster_number']
                num_in_clust = row['number_tcrs_in_cluster']
                mult = row['total_multiplicity']
                is_assoc = row['is_cancer_assoc']
                if clust_num in temp_dict.keys():
                    temp_dict[clust_num]['contrib_pts'] += 1
                    temp_dict[clust_num]['num_in_clust'] += num_in_clust
                    temp_dict[clust_num]['mult'] += mult
                else:
                    temp_dict[clust_num] = {'contrib_pts': 1, 
                                          'num_in_clust': num_in_clust, 
                                          'mult': mult, 'is_assoc': is_assoc}
        df_dict = {'cluster_number': [], 'contributing_patients': [], 
                   'num_tcrs_in_clust': [], 'total_reads': [], 
                   'is_cancer_assoc': []}
        for k, v in temp_dict.items():
            df_dict['cluster_number'].append(k)
            df_dict['contributing_patients'].append(v['contrib_pts'])
            df_dict['num_tcrs_in_clust'].append(v['num_in_clust'])
            df_dict['total_reads'].append(v['mult'])
            df_dict['is_cancer_assoc'].append(v['is_assoc'])
        output_df = pd.DataFrame(df_dict)
        output_df.to_csv(
            join(group_assoc_summaries_dir, labels_reformat[i] + '.csv'), 
            index=False, sep='\t')
            

def assoc_btwn_groups(group_assoc_summaries_dir, btwn_group_summary_dir):
    summaries = sorted(listdir(group_assoc_summaries_dir))
    summ_names = [s[:-4] for s in summaries]

    long_int_df = pd.read_csv(join(group_assoc_summaries_dir, summaries[0]), 
                           sep='\t', low_memory=False)
    no_nact_df = pd.read_csv(join(group_assoc_summaries_dir, summaries[1]), 
                          sep='\t', low_memory=False)
    short_int_df = pd.read_csv(join(group_assoc_summaries_dir, summaries[2]), 
                            sep='\t', low_memory=False)
    df_arr = [long_int_df, no_nact_df, short_int_df]
    
    idx_combos = [[0, 1], [1, 2], [0, 2]]

    for combo in idx_combos:
        idx1 = combo[0]
        idx2 = combo[1]
        df1 = df_arr[idx1]
        df2 = df_arr[idx2]
        summ_name_1 = summ_names[idx1]
        summ_name_2 = summ_names[idx2]

        rename1 = {'contributing_patients':  f'contributing_patients_{summ_name_1}', 
                    'num_tcrs_in_clust': f'num_tcrs_in_clust_{summ_name_1}',
                    'total_reads': f'total_reads_{summ_name_1}'}
        df1 = df1.rename(columns=rename1)

        rename2 = {'contributing_patients':  f'contributing_patients_{summ_name_2}', 
                    'num_tcrs_in_clust': f'num_tcrs_in_clust_{summ_name_2}',
                    'total_reads': f'total_reads_{summ_name_2}'}
        df2 = df2.rename(columns=rename2)
        output_df = pd.merge(df1.iloc[:,:-1], df2, on='cluster_number')
        
        fname = summ_name_1 + '_vs_' + summ_name_2 + '.csv'
        full_path = join(btwn_group_summary_dir, fname)
        output_df.to_csv(full_path, index=False, sep='\t')

        # set1 = set(df1['cluster_number'])
        # set2 = set(df2['cluster_number'])
        # intsect = set1.intersection(set2)

        # filter_fn = lambda x: x['cluster_number'] in intsect

        # slice1 = df1[df1.apply(filter_fn, axis=1)]
        # sorted1 = slice1.sort_values('cluster_number')
        # rename1 = {'contributing_patients':  f'contributing_patients_{summ_name_1}', 
        #             'num_tcrs_in_clust': f'num_tcrs_in_clust_{summ_name_1}',
        #             'total_reads': f'total_reads_{summ_name_1}', 
        #             'is_cancer_assoc': f'is_cancer_assoc_{summ_name_1}'}
        # formatted1 = sorted1.rename(columns=rename1)
        # formatted1.reset_index(inplace=True, drop=True)
        
        # slice2 = df2[df2.apply(filter_fn, axis=1)]
        # sorted2 = slice2.sort_values('cluster_number')
        # rename2 = {'contributing_patients':  f'contributing_patients_{summ_name_2}', 
        #             'num_tcrs_in_clust': f'num_tcrs_in_clust_{summ_name_2}',
        #             'total_reads': f'total_reads_{summ_name_2}', 
        #             'is_cancer_assoc': f'is_cancer_assoc_{summ_name_2}'}
        # formatted2 = sorted2.iloc[:, 1:].rename(columns=rename2)
        # formatted2.reset_index(inplace=True, drop=True)
        # full_df = pd.concat([formatted1, formatted2], axis=1)


    # # Long vs no nact
    # long_set = set(long_int_df['cluster_number'])
    # no_nact_set = set(no_nact_df['cluster_number'])
    # long_no_intsect = long_set.intersection(no_nact_set)
    
    # filter_fn = lambda x: x['cluster_number'] in long_no_intsect
    
    # long_slice = long_int_df[long_int_df.apply(filter_fn, axis=1)]
    # long_sorted = long_slice.sort_values('cluster_number')
    # long_rename = {'contributing_patients':  'contributing_patients_long_nact', 
    #                'num_tcrs_in_clust': 'num_tcrs_in_clust_long_nact',
    #                'total_reads': 'total_reads_long_nact', 
    #                'is_cancer_assoc': 'is_cancer_assoc_long_nact'}
    # long_formatted = long_sorted.iloc[:, 1:].rename(columns=long_rename)
    # long_formatted.reset_index(inplace=True, drop=True)
    
    # no_slice = no_nact_df[no_nact_df.apply(filter_fn, axis=1)]
    # no_sorted = no_slice.sort_values('cluster_number')
    # no_rename = {'contributing_patients':  'contributing_patients_no_nact', 
    #               'num_tcrs_in_clust': 'num_tcrs_in_clust_no_nact',
    #               'total_reads': 'total_reads_no_nact', 
    #               'is_cancer_assoc': 'is_cancer_assoc_no_nact'}
    # no_formatted = no_sorted.iloc[:, 1:].rename(columns=no_rename)
    # no_formatted.reset_index(inplace=True, drop=True)
    # pd.concat([long_formatted, no_formatted], axis=1)
    