"""
Utility function for tests.
"""

def assertAlmostEqualList(testcase, actualList, expectedList, delta):
    assert len(actualList) == len(expectedList)
    for actual, expected in zip(actualList, expectedList):
        testcase.assertAlmostEqual(actual, expected, delta=delta,
                                   msg=str(actualList) + " != " + str(expectedList))

def get_first_index_containing_substring(string_list, substring, default_index=-1):  
    """  
    Returns the index of the first item in a list of strings that contains the given substring.  
    This search is case sensitive.   
    :param string_list: The list of strings to search.  
    :param substring: The substring to match.  
    :param default_index: Index to return if substring is not found.  
    :return: The index of the first item containing the substring, or default_index if not found.  
    """  
    for index, s in enumerate(string_list):  
        if substring in s:  
            return index  
    return default_index  

def check_annotation_term_ids(term_ids: list):  
    """  
    Check that annotation terms have "UBERON", "ILX" or "" (empty string) in the primary term id, and  
    that any supplied term ids are in the order UBERON < ILX < FMA.  
    String comparisons are case senstitive.  
    :param term_ids: List of string annotation names or term identifiers.  
    :return: True if annotation term ids are valid and in the required order, otherwise False.  
    """  
    first_id = term_ids[1]  
    if ("UBERON" not in first_id) and ("ILX" not in first_id) and (first_id != ""): 
        return False  
    term_ids_count = len(term_ids)  
    uberon_index = get_first_index_containing_substring(term_ids, "UBERON", -1)  
    ilx_index = get_first_index_containing_substring(term_ids, "ILX", term_ids_count)  
    fma_index = get_first_index_containing_substring(term_ids, "FMA", 2 * term_ids_count)  
    if (uberon_index >= 0) and ((ilx_index < uberon_index) or (fma_index < uberon_index)):  
        return False  
    if (ilx_index < term_ids_count) and (fma_index < ilx_index):  
        return False  
    return True  


