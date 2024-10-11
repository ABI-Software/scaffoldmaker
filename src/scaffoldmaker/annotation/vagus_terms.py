"""
Common resource for vagus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
vagus_marker_terms = [
    # anatomical landmarks
    ("level of superior border of jugular foramen on the vagus nerve", "ILX:0794617"),
    ("right level of superior border of jugular foramen on the vagus nerve", "ILX:0794618"),
    ("left level of superior border of jugular foramen on the vagus nerve", "ILX:0794619"),
    ("level of inferior border of jugular foramen on the vagus nerve", "ILX:0794620"),
    ("right level of inferior border of jugular foramen on the vagus nerve", "ILX:0794621"),
    ("left level of inferior border of jugular foramen on the vagus nerve", "ILX:0794622"),
    ("level of inferior border of cranium on the vagus nerve", "ILX:0794623"),
    ("right level of inferior border of cranium on the vagus nerve", "ILX:0794624"),
    ("left level of inferior border of cranium on the vagus nerve", "ILX:0794625"),
    ("level of C1 transverse process on the vagus nerve", "ILX:0794626"),
    ("right level of C1 transverse process on the vagus nerve", "ILX:0794627"),
    ("left level of C1 transverse process on the vagus nerve", "ILX:0794628"),
    ("level of greater horn of hyoid on the vagus nerve", "ILX:0794629"),
    ("right level of greater horn of hyoid on the vagus nerve", "ILX:0794630"),
    ("left level of greater horn of hyoid on the vagus nerve", "ILX:0794631"),
    ("level of laryngeal prominence on the vagus nerve", "ILX:0794632"),
    ("right level of laryngeal prominence on the vagus nerve", "ILX:0794633"),
    ("left level of laryngeal prominence on the vagus nerve", "ILX:0794634"),
    ("level of angle of the mandible on the vagus nerve", "ILX:0794635"),
    ("right level of angle of the mandible on the vagus nerve", "ILX:0794636"),
    ("left level of angle of the mandible on the vagus nerve", "ILX:0794637"),
    ("level of carotid bifurcation on the vagus nerve", "ILX:0794638"),
    ("right level of carotid bifurcation on the vagus nerve", "ILX:0794639"),
    ("left level of carotid bifurcation on the vagus nerve", "ILX:0794640"),
    ("level of superior border of the clavicle on the vagus nerve", "ILX:0794641"),
    ("right level of superior border of the clavicle on the vagus nerve", "ILX:0794642"),
    ("left level of superior border of the clavicle on the vagus nerve", "ILX:0794643"),
    ("level of jugular notch on the vagus nerve", "ILX:0794644"),
    ("right level of jugular notch on the vagus nerve", "ILX:0794645"),
    ("left level of jugular notch on the vagus nerve", "ILX:0794646"),
    #("level of carina", "None"),
    #("right level of carina", "None"),
    #("left level of carina", "None"),
    ("level of sternal angle on the vagus nerve", "ILX:0794647"),
    ("right level of sternal angle on the vagus nerve", "ILX:0794648"),
    ("left level of sternal angle on the vagus nerve", "ILX:0794649"),
    ("level of 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794650"),  # !!!
    ("right level of 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794651"),  # !!!
    ("left level of 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794652"),  # !!!
    ("level of esophageal hiatus on the vagus nerve", "ILX:0794653"),
    ("right level of esophageal hiatus on the vagus nerve", "ILX:0794654"),
    ("left level of esophageal hiatus on the vagus nerve", "ILX:0794655"),
    ("level of aortic hiatus on the vagus nerve", "ILX:0794656"),
    ("right level of aortic hiatus on the vagus nerve", "ILX:0794657"),
    ("left level of aortic hiatus on the vagus nerve", "ILX:0794658"),
]

vagus_branch_terms = [
    # branch names, based on lit.review
    # branch names with no indication of the side
    ("meningeal branch of vagus nerve", "FMA:6231"),
    ("communicating branch of vagus nerve with glossopharyngeal nerve", "FMA:6233"),
    ("auricular branch of vagus nerve", "FMA:6232"),
    ("pharyngeal branch of vagus nerve", "UBERON:0000929"),
    ("lingual branch of vagus nerve", "FMA:6235"),
    ("branch of vagus nerve to carotid body", "FMA:6237"),
    ("communicating branch of superior cervical ganglion with vagus nerve", "FMA:6901"),
    ("superior laryngeal nerve", "UBERON:0011326"),
    ("internal laryngeal nerve", "FMA:6240"),
    ("external branch of superior laryngeal nerve", "FMA:6243"),
    ("communicating branch of internal laryngeal nerve with recurrent laryngeal nerve", "FMA:53544"),
    ("communicating branch of external laryngeal nerve with superior cardiac nerve", "FMA:6708"),
    ("pulmonary branch of vagus nerve", "FMA:65515"),
    ("bronchial branch of vagus nerve", "FMA:6247"),
    ("superior cervical cardiac branch of vagus nerve", "FMA:6244"),
    ("recurrent laryngeal nerve", "UBERON:0003716"),
    ("extra-laryngeal branch of recurrent laryngeal nerve to larynx", "FMA:6710"),
    ("oesophageal branch of recurrent laryngeal nerve", "FMA:6248"),
    ("tracheal branch of recurrent laryngeal nerve", "FMA:6249"),
    ("inferior laryngeal nerve", "UBERON:0003716"), # same code as for recurrent laryngeal nerve
    ("communicating branch of recurrent laryngeal nerve with internal laryngeal nerve", "FMA:53526"),
    ("inferior cervical cardiac branch of vagus nerve", "FMA:75530"),
    ("thoracic cardiac branch of vagus nerve", "FMA:53601"),

    # right vagus branches
    ("right vagus nerve", "FMA:6219"),
    ("right vagus X nerve trunk", "UBERON:0035021"),
    ("meningeal branch of right vagus nerve", "FMA:53541"),
    ("recurrent meningeal branch of right vagus nerve", "None"),
    ("communicating branch of right vagus nerve with right glossopharyngeal nerve", "FMA:52559"),
    ("auricular branch of right vagus nerve", "FMA:53534"),
    ("communicating branch of auricular branch of right vagus nerve with right facial nerve", "FMA:53587"),
    ("communicating branch of auricular branch of right vagus nerve with right posterior auricular nerve", "FMA:53589"),
    ("pharyngeal branch of right vagus nerve to pharyngeal nerve plexus", "FMA:53635"),
    ("lingual branch of right vagus nerve", "FMA:53633"),
    ("communicating branch of pharyngeal branch of right vagus nerve with superior cervical ganglion", "None"),
    ("branch of right vagus nerve to carotid body", "FMA:53606"),
    ("communicating branch of superior cervical ganglion with right vagus nerve", "ILX:0794055"),
    ("right superior laryngeal nerve", "FMA:53530"),
    ("right internal laryngeal nerve", "FMA:53539"),
    ("right external laryngeal nerve", "FMA:53537"),
    ("superior branch of right internal laryngeal nerve", "FMA:53575"), # generic term???
    ("middle branch of right internal laryngeal nerve", "None"),
    ("inferior branch of right internal laryngeal nerve", "FMA:53581"), # generic term???
    ("communicating branch of right internal laryngeal nerve with right recurrent laryngeal nerve", "FMA:53571"),
    ("communicating branch of right external laryngeal nerve with right superior cardiac nerve", "FMA:53561"),
    ("cardiopulmonary branch A of right vagus nerve", "ILX:0794154"),
    ("cardiopulmonary branch B of right vagus nerve", "ILX:0794155"),
    ("pulmonary branch of right vagus nerve", "FMA:6671"),
    ("pulmonary branch A of right vagus nerve", "None"),
    ("pulmonary branch B of right vagus nerve", "None"),
    ("pulmonary branch C of right vagus nerve", "None"),
    ("pulmonary branch D of right vagus nerve", "None"),
    ("pulmonary branch E of right vagus nerve", "None"),
    # === temporary branches used in the japanese dataset
    ("right A lateral pulmonary branch of vagus nerve", "ILX:0794278"),
    ("right B lateral pulmonary branch of vagus nerve", "ILX:0794279"),
    ("posterior A gastric branch of vagus nerve", "ILX:0794398"),
    ("posterior B gastric branch of vagus nerve", "ILX:0794399"),
    ("posterior C gastric branch of vagus nerve", "ILX:0794400"),
    # ===
    ("bronchial branch of right vagus nerve", "FMA:53613"),
    ("superior cervical cardiac branch of right vagus nerve", "FMA:53598"),
    ("right recurrent laryngeal nerve", "UBERON:0011767"),
    ("extra-laryngeal branch of right recurrent laryngeal nerve to larynx", "None"),
    ("branch of right recurrent laryngeal nerve to muscle of larynx", "None"),
    ("oesophageal branch of right recurrent laryngeal nerve", "FMA:53608"),
    ("tracheal branch of right recurrent laryngeal nerve", "FMA:53610"),
    ("right inferior laryngeal nerve", "UBERON:0011767"), # same code as for left recurrent laryngeal nerve
    ("external branch of right inferior laryngeal nerve", "None"),
    ("anterior branch of right inferior laryngeal nerve", "None"),
    ("posterior branch of right inferior laryngeal nerve", "None"),
    ("communicating branch of right recurrent laryngeal nerve with superior cervical ganglion", "None"),
    ("communicating branch of right recurrent laryngeal nerve with right internal laryngeal nerve", "FMA:53528"),
    ("inferior cervical cardiac branch of right recurrent laryngeal nerve", "ILX:0794235"),
    ("inferior cervical cardiac branch of right vagus nerve", "FMA:6713"),
    ("thoracic cardiac branch of right vagus nerve", "FMA:53604"),
    ("right cardiac branch to deep cardiac nerve plexus", "FMA:6711"),
    ("branch of right vagus nerve to oesophageal nerve plexus", "ILX:0794299"),

    # posterior vagus, a continuation of right vagus
    ("Posterior esophageal vagus trunk", "ILX:0794858"),
    ("Celiac branch of posterior vagal trunk", "FMA:6667"),
    ("Greater posterior gastric nerve", "FMA:6689"),
    ("Pyloric branch of greater posterior gastric nerve", "FMA:6677"),

    # left vagus branches
    ("left vagus nerve", "FMA:6220"),
    ("left vagus X nerve trunk", "UBERON:0035020"),
    ("meningeal branch of left vagus nerve", "FMA:53542"),
    ("recurrent meningeal branch of left vagus nerve", "None"),
    ("communicating branch of left vagus nerve with left glossopharyngeal nerve", "FMA:53560"),
    ("auricular branch of left vagus nerve", "FMA:53535"),
    ("communicating branch of auricular branch of left vagus nerve with left facial nerve", "FMA:53588"),
    ("communicating branch of auricular branch of left vagus nerve with left posterior auricular nerve", "FMA:53590"),
    ("pharyngeal branch of left vagus nerve to pharyngeal nerve plexus", "FMA:53636"),
    ("lingual branch of left vagus nerve", "FMA:53634"),
    ("communicating branch of pharyngeal branch of left vagus nerve with superior cervical ganglion", "None"),
    ("branch of left vagus nerve to carotid body", "FMA:53607"),
    ("communicating branch of superior cervical ganglion with left vagus nerve", "ILX:0794061"),
    ("left superior laryngeal nerve", "FMA:53536"),
    ("left internal laryngeal nerve", "FMA:53540"),
    ("left external laryngeal nerve", "FMA:53538"),
    ("superior branch of left internal laryngeal nerve", "FMA:53576"), # generic term???
    ("middle branch of left internal laryngeal nerve", "None"),
    ("inferior branch of left internal laryngeal nerve", "FMA:53582"), # generic term???
    ("communicating branch of left internal laryngeal nerve with left recurrent laryngeal nerve", "FMA:53572"),
    ("communicating branch of left external laryngeal nerve with left superior cardiac nerve", "FMA:53562"),
    ("cardiopulmonary branch A of left vagus nerve", "ILX:0794160"),
    ("cardiopulmonary branch B of left vagus nerve", "ILX:0794161"),
    ("pulmonary branch of left vagus nerve", "FMA:6679"),
    ("pulmonary branch A of left vagus nerve", "None"),
    ("pulmonary branch B of left vagus nerve", "None"),
    ("pulmonary branch C of left vagus nerve", "None"),
    ("pulmonary branch D of left vagus nerve", "None"),
    # === temporary branches used in the japanese dataset
    ("left A medial pulmonary branch of vagus nerve", "ILX:0794266"),
    ("left B medial pulmonary branch of vagus nerve", "ILX:0794267"),
    # ===
    ("bronchial branch of left vagus nerve", "FMA:53558"),
    ("superior cervical cardiac branch of left vagus nerve", "FMA:53599"),
    ("left recurrent laryngeal nerve", "UBERON:0011766"),
    ("extra-laryngeal branch of left recurrent laryngeal nerve to larynx", "None"),
    ("branch of left recurrent laryngeal nerve to muscle of larynx", "None"),
    ("oesophageal branch of left recurrent laryngeal nerve", "FMA:53609"),
    ("tracheal branch of left recurrent laryngeal nerve", "FMA:53611"),
    ("left inferior laryngeal nerve", "UBERON:0011766"),
    ("external branch of left inferior laryngeal nerve", "None"),
    ("anterior branch of left inferior laryngeal nerve", "None"),
    ("posterior branch of left inferior laryngeal nerve", "None"),
    ("communicating branch of left recurrent laryngeal nerve with superior cervical ganglion", "None"),
    ("communicating branch of left recurrent laryngeal nerve with left internal laryngeal nerve", "FMA:53529"),
    ("inferior cervical cardiac branch of left recurrent laryngeal nerve", "ILX:0794244"),
    ("inferior cervical cardiac branch of left vagus nerve", "FMA:6714"),
    ("thoracic cardiac branch of left vagus nerve", "FMA:53605"),
    ("left cardiac branch to deep cardiac nerve plexus", "None"),
    ("branch of left vagus nerve to oesophageal nerve plexus", "ILX:0794310"),

    # anterior vagus, a continuation of left vagus
    ("anterior esophageal vagus trunk", "ILX:0794854"),
    ("hepatic branch of anterior vagal trunk", "FMA:6666"),
    ("greater anterior gastric nerve", "FMA:6684"),
    ("branch of greater anterior gastric nerve to coeliac nerve plexus", "FMA:53675"),

]

def access_vagus_marker_terms():
    return vagus_marker_terms

def access_vagus_branch_terms():
    return vagus_branch_terms

# def get_vagus_branch_term(name : str):
#     """
#     Find term by matching name to any identifier held for a term.
#     Raise exception if name not found.
#     :return ( preferred name, preferred id )
#     """
#     for term in vagus_branch_terms:
#         if name in term:
#             return ( term[0], term[1] )
#     raise NameError("Vagus annotation term '" + name + "' not found.")

def get_vagus_marker_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in vagus_marker_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Vagus annotation term '" + name + "' not found.")

def marker_name_in_terms(name : str):
    for term in vagus_marker_terms:
        if name in term:
            return True
    return False

def get_vagus_branch_term(name, vagus_terms):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in vagus_terms:
        if name in term:
            return ( term[0], term[1] )
    return ( name, "")
