import pandas as pd
import numpy as np

def parse_sample_overview(samp_overview_path, saved_df_path):
    ov_df = pd.read_csv(samp_overview_path, sep='\t')
    ov_df = ov_df.sort_values(by=['sample_name'])

    # Splits tag column into list of strings
    tag_split_lambda = lambda x: x.lower().split(', ')
    ov_df['sample_tags'] = ov_df['sample_tags'].apply(tag_split_lambda)

    # Pulls out relevant values from list of strings in sample tag column
    ov_df['nact_intvl'] = ov_df['sample_tags'].apply(_extract_intvl_nact)
    ov_df['survival_months'] = ov_df['sample_tags'].apply(_extract_survival)
    ov_df['age'] = ov_df['sample_tags'].apply(_extract_age)
    ov_df['sex'] = ov_df['sample_tags'].apply(_extract_sex)
    ov_df['biopsy_location'] = ov_df['sample_tags'].apply(_extract_biop_loc)
    ov_df['interval_group_num'] = ov_df['sample_tags'].apply(_extract_int_grp)
    ov_df['rel_prod_rearr'] = \
        ov_df['productive_rearrangements'] / ov_df['productive_templates']
    ov_df['rel_total_rearr'] = \
        ov_df['total_rearrangements'] / ov_df['total_templates']
    ov_df['frac_prod_rearr'] = \
        ov_df['productive_rearrangements'] / ov_df['total_rearrangements']
    group_label_list = ['No NACT', 'Short Interval', 'Long Interval']
    ov_df['group_label'] = \
        ov_df['interval_group_num'].apply(lambda x: group_label_list[x])

    ov_df.to_csv(saved_df_path, sep="\t", index=False)
    

def _containing_substr(str_list, substr):
    """
    Given a list of strings, returns the first string in the list that contains
    a given substring.

    Args:
        str_list (list(str)): list of strings to search
        substr (str): substring to search for
    
    Returns:
        str: first string containing substring
    """
    for str_i in str_list:
        if substr in str_i:
            return str_i


def _extract_intvl_nact(list_str):
    """
    Given a list of strings, pull out the interval value of NACT treatment.
    Listed as np.NaN for patients who didn't recieve treatment.

    Args:
        list_str (list(str)): list of strings parsed from tag column of
        original datafile

    Returns:
        float: interval of NACT treatment
    """
    substr = _containing_substr(list_str, ' week nact interval')
    extracted_num = substr.split(' ')[0]
    if extracted_num == 'na':
        return np.NaN
    else:
        return float(extracted_num)


def _extract_survival(list_str):
    """
    Given a list of strings, pull out the number of months of survival.

    Args:
        list_str (list(str)): list of strings parsed from tag column of
        original datafile

    Returns:
        int: number of months of survival
    """
    substr = _containing_substr(list_str, ' month survival')
    extracted_num = substr.split(' ')[0]
    return int(extracted_num)


def _extract_age(list_str):
    """
    Given a list of strings, pull out the patient age.

    Args:
        list_str (list(str)): list of strings parsed from tag column of
        original datafile

    Returns:
        int: age of patient
    """
    substr = _containing_substr(list_str, ' years')
    extracted_num = substr.split(' ')[0]
    return int(extracted_num)


def _extract_sex(list_str):
    """
    Given a list of strings, pull out the patient sex.

    Args:
        list_str (list(str)): list of strings parsed from tag column of
        original datafile

    Returns:
        str: patient sex
    """
    return _containing_substr(list_str, 'male')


def _extract_biop_loc(list_str):
    """
    Given a list of strings, pull out the biopsy location.

    Args:
        list_str (list(str)): list of strings parsed from tag column of
        original datafile

    Returns:
        str: biopsy location
    """
    str_avoid_list = [' week nact interval', ' month survival', 
                      ' years', 'male', 'nact interval group ']
    for str_i in list_str:
        if not any(x in str_i for x in str_avoid_list):
            return str_i


def _extract_int_grp(list_str):
    """
    Given a list of strings, pull out the NACT treatment interval group number.

    Args:
        list_str (list(str)): list of strings parsed from tag column of
        original datafile

    Returns:
        int: NACT treatment group number
    """
    substr = _containing_substr(list_str, 'nact interval group ')
    extracted_num = substr.split(' ')[3]
    return int(extracted_num)