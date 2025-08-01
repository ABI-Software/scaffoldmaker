"""
Utilities for working with annotations.
"""
import logging


logger = logging.getLogger(__name__)


def annotation_term_id_to_url(term):
    """
    Convert short forms of annotation term ids into full urls e.g.:
    "UBERON:0000948" --> "http://purl.obolibrary.org/obo/UBERON_0000948"
    "ILX:0732254" --> "http://uri.interlex.org/base/ilx_0732254"
    "FMA:7088" --> "http://purl.org/sig/ont/fma/fma7088"
    Other namespaces should be added as needed, but until then are returned as is with a warning.
    :param term_name:
    :param term_id: String short form of annotation term id NAMESPACE:#.
    :return: term with id converted to url idfpossible.
    """
    if not term[1]:
        logger.warning("annotation_term_id_to_url:  Missing term id for annotation '" + term[0] + "'")
        return term
    if ":" not in term[1]:
        logger.error("annotation_term_id_to_url:  Invalid annotation term id ('" + term[0] + "', '" + term[1] + "')")
        return term
    raw_prefix, number_str = term[1].split(":")
    prefix = raw_prefix.upper()
    if prefix == "UBERON":
        url_prefix = "http://purl.obolibrary.org/obo/UBERON_"
    elif prefix == "ILX":
        url_prefix = "http://uri.interlex.org/base/ilx_"
    elif prefix == "FMA":
        url_prefix = "http://purl.org/sig/ont/fma/fma"
    else:
        logger.warning("annotation_term_id_to_url:  No url known for term ('" + term[0] + "', '" + term[1] + "')")
        return term
    return (term[0], url_prefix + number_str)
