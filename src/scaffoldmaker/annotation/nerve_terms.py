"""
Common resource for nerve centreline annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
nerve_terms = [
    ('vagus nerve', 'UBERON:0001759', 'ILX:0733375'),
    ("trigeminal nucleus", "UBERON:0002925"),
    ("nucleus of solitary tract", "UBERON:0009050"),
    ("dorsal motor nucleus of vagus nerve", "UBERON:0002870"),
    ("nucleus ambiguus", "UBERON:0001719"),
    ('point_1', ''),
    ('ardell_4_branching_point', ''),
    ("superior vagus X ganglion", "UBERON:0005364"),
    ('bolser_X_1', ''),
    ('ardell_3_branching_point', ''),
    ('point_3', ''),
    ("dura mater in the posterior cranial fossa", ""),
    ("inferior vagus X ganglion", "UBERON:0005363"),
    ('point_6', ''),
    ('point_7', ''),
    ('point_8', ''),
    ('point_9', ''),
    ('point_10', ''),
    ("tympanic membrane", "UBERON:0002364"),
    ("floor of external acoustic meatus", ""),
    ("posterior wall of external acoustic meatus", ""),
    ('point_11', ''),
    ('point_13', ''),
    ("pharynx)", "UBERON:0006562"),
    ("back part of the tongue*", "UBERON:0010033"),
    ('point_14', ''),
    ("carotid body", "UBERON:0001629"),
    ('point_16', ''),
    ('point_18', ''),
    ('point_19', ''),
    ('point_20', ''),
    ("mucosa of pharynx", "FMA:55031"),
    ("epiglottic vallecula", "UBERON:0013165"),
    ("epiglottis", "UBERON:0000388"),
    ("Wall of larynx", "ILX:0738326"),
    ('point_21', ''),
    ("aryepiglottic fold", "UBERON:0014385"),
    ("arytenoideus", ""),
    ("mucosa of arytenoid cartilage", ""),
    ('point_37', ''),
    ('bolser_X-5', ''),
    ('bolser_X-4', ''),
    ('external_laryngeal_n_branching_point', ''),
    ("pharyngeal nerve plexus", "UBERON:0011325"),
    ('bolser_X-14', ''),
    ("inferior pharyngeal constrictor", "UBERON:0001570"),
    ("cricothyroid muscle", "UBERON:0001566"),
    ('bolser_X-15', ''),
    ('point_41', ''),
    ('point_22', ''),
    ('point_24', ''),
    ('point_25', ''),
    ("cardiac nerve plexus", "UBERON:0002008"),
    ('point_26', ''),
    ('point_28', ''),
    ("mucosa of larynx", "UBERON:0001824"),
    ('point_29', ''),
    ("inferior pharyngeal constrictor", "UBERON:0001570"),
    ("ceratocricoid", ""),
    ("lateral crico-arytenoid", "UBERON:0008573"),
    ("oblique arytenoid", "UBERON:0008575"),
    ("posterior crico-arytenoid", "UBERON:0008572"),
    ("hyro-arytenoid", "FMA:46588"),
    ("transverse arytenoid", "UBERON:0008574"),
    ("vocalis muscle", "UBERON:0008577"),
    ("vocal_cords", ""),
    ("Laryngeal mechanoreceptors", "ILX:0738317"),
    ("epiglottis", "UBERON:0000388"),
    ('point_38', ''),
    ('bolser_X-2', ''),
    ('bolser_X-3', ''),
    ('point_30', ''),
    ("pulmonary nerve plexus", "UBERON:0002009"),
    ('point_31', ''),
    ("bronchus smooth muscle", "UBERON:0004242"),
    ('point_32', ''),
    ('point_33', ''),
    ("esophageal nerve plexus", "FMA:6225"),
    ("esophagus", "UBERON:0001043"),
    ('point_34', ''),
    ("trigeminal nucleus", "UBERON:0002925"),
    ("pancreas", "UBERON:0001264"),
    ('point_35', ''),
    ("left gastric nerve plexus", "FMA:6633"),
    ('point_36', ''),
    ("kidney", "UBERON:0002113"),
    # flatmap centrelines
    ('ardell_3', ''),
    ('ardell_4', ''),
    ('bolser_1', ''),
    ('bolser_2', ''),
    ('bolser_3', ''),
    ('bolser_4', ''),
    ('bolser_5', ''),
    ('bolser_14', ''),
    ('bolser_15', ''),
    ('n_1', ''),
    ('n_2', ''),
    ('n_3', ''),
    ('n_4', ''),
    ('n_5', ''),
    ('n_6', ''),
    ('n_68', ''),
    ('n_7', ''),
    ('n_8', ''),
    ('n_9', ''),
    ('n_10', ''),
    ('n_11', ''),
    ('n_12', ''),
    ('n_13', ''),
    ('n_14', ''),
    ('n_15', ''),
    ('n_16', ''),
    ('n_17', ''),
    ('n_18', ''),
    ('n_19', ''),
    ('n_20', ''),
    ('n_21', ''),
    ('n_22', ''),
    ('n_23', ''),
    ('n_24', ''),
    ('n_25', ''),
    ('n_26', ''),
    ('n_27', ''),
    ('n_28', ''),
    ('n_29', ''),
    ('n_30', ''),
    ('n_31', ''),
    ('n_32', ''),
    ('n_33', ''),
    ('n_34', ''),
    ('n_35', ''),
    ('n_36', ''),
    ('n_37', ''),
    ('n_38', ''),
    ('n_39', ''),
    ('n_40', ''),
    ('n_41', ''),
    ('n_42', ''),
    ('n_43', ''),
    ('n_44', ''),
    ('n_45', ''),
    ('n_46', ''),
    ('n_47', ''),
    ('n_48', ''),
    ('n_49', ''),
    ('n_50', ''),
    ('n_51', ''),
    ('n_52', ''),
    ('n_53', ''),
    ('n_54', ''),
    ('n_55', ''),
    ('n_56', ''),
    ('n_57', ''),
    ('n_58', ''),
    ('n_59', ''),
    ('n_60', ''),
    ('n_61', ''),
    ('n_62', ''),
    ('n_63', ''),
    ('n_64', ''),
    ('n_65', ''),
    ('n_66', ''),
    ('n_67', ''),
    ]

def get_nerve_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in nerve_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Nerve annotation term '" + name + "' not found.")
