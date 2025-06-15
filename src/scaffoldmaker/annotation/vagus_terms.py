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
    # ("level of carina", "None"),
    # ("right level of carina", "None"),
    # ("left level of carina", "None"),
    ("level of sternal angle on the vagus nerve", "ILX:0794647"),
    ("right level of sternal angle on the vagus nerve", "ILX:0794648"),
    ("left level of sternal angle on the vagus nerve", "ILX:0794649"),
    ("level of 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794650"),
    ("right level of 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794651"),
    ("left level of 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794652"),
    ("level of esophageal hiatus on the vagus nerve", "ILX:0794653"),
    ("right level of esophageal hiatus on the vagus nerve", "ILX:0794654"),
    ("left level of esophageal hiatus on the vagus nerve", "ILX:0794655"),
    ("level of aortic hiatus on the vagus nerve", "ILX:0794656"),
    ("right level of aortic hiatus on the vagus nerve", "ILX:0794657"),
    ("left level of aortic hiatus on the vagus nerve", "ILX:0794658"),
    # vagus built-in annotations
    ("vagus centroid", ""),
    ("vagus epineureum", ""),
    ("vagus anterior line", "")
]

vagus_branch_terms = [
    # branch names, based on lit.review
    # branch names with no indication of the side
    ("cervical trunk", "ILX:0789914"),
    ("thoracic trunk", "ILX:0784729"),
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
    ("esophageal branch of recurrent laryngeal nerve", "FMA:6248"),
    ("tracheal branch of recurrent laryngeal nerve", "FMA:6249"),
    ("inferior laryngeal nerve", "UBERON:0003716"),  # same code as for recurrent laryngeal nerve
    ("communicating branch of recurrent laryngeal nerve with internal laryngeal nerve", "FMA:53526"),
    ("inferior cervical cardiac branch of vagus nerve", "FMA:75530"),
    ("thoracic cardiac branch of vagus nerve", "FMA:53601"),
    ("esophageal trunk", "ILX:0794853"),

    # right vagus branches
    ("right vagus nerve", "FMA:6219", "ILX:0789705"),
    ("right vagus X nerve trunk", "UBERON:0035021", "ILX:0730515"),
    ("right cervical trunk", "ILX:0794141"),
    ("right thoracic trunk", "ILX:0786664"),
    ("right meningeal branch of right vagus nerve", "FMA:53541", "ILX:0785804"),
    ("right branch between vagus nerve and glossopharyngeal nerve", "FMA:53559", "ILX:0790506"),
    ("right auricular branch of right vagus nerve", "FMA:53534", "ILX:0785879"),
    ("communicating branch of auricular branch of right vagus nerve with right facial nerve", "FMA:53587", "ILX:0791102"),
    ("communicating branch of auricular branch of right vagus nerve with right posterior auricular nerve", "FMA:53589", "ILX:0786954"),
    ("right pharyngeal branch of right vagus nerve to pharyngeal nerve plexus", "FMA:53635", "ILX:0792814"),
    ("lingual branch of right vagus nerve", "FMA:53633", "ILX:0787083"),
    ("right pharyngeal branch of right vagus nerve to superior cervical ganglion", "ILX:0795066"),
    ("branch of right vagus nerve to carotid body", "FMA:53606", "ILX:0790433"),
    ("right branch between vagus nerve and superior cervical ganglion", "ILX:0794055"),
    ("right superior laryngeal nerve", "FMA:53530", "ILX:0787738"),
    ("right internal laryngeal nerve", "FMA:53539", "ILX:0788164"),
    ("right external laryngeal nerve", "FMA:53537", "ILX:0792879"),
    ("superior branch of right internal laryngeal nerve", "FMA:53575", "ILX:0784524"),
    ("middle branch of right internal laryngeal nerve", "None"),  # NA pending
    ("inferior branch of right internal laryngeal nerve", "FMA:53581", "ILX:0789165"),
    ("communicating branch of right internal laryngeal nerve with right recurrent laryngeal nerve", "FMA:53571", "ILX:0791006"),
    ("communicating branch of right external laryngeal nerve with right superior cardiac nerve", "FMA:53561", "ILX:0787405"),
    ("right A cervical cardiopulmonary branch of vagus nerve", "ILX:0794154"),
    ("right B cervical cardiopulmonary branch of vagus nerve", "ILX:0794155"),
    ("right pulmonary branch of vagus nerve", "FMA:6671", "ILX:0787735"),
    ("right pulmonary branch A of the vagus nerve", "ILX:0795074"),
    ("right pulmonary branch B of the vagus nerve", "ILX:0795075"),
    ("right pulmonary branch C of the vagus nerve", "ILX:0795076"),
    ("right pulmonary branch D of the vagus nerve", "ILX:0795077"),
    ("right pulmonary branch E of the vagus nerve", "ILX:0795078"),
    # === temporary branches used in the japanese dataset
    ("right A lateral pulmonary branch of vagus nerve", "ILX:0794278"),
    ("right B lateral pulmonary branch of vagus nerve", "ILX:0794279"),
    ("posterior A gastric branch of vagus nerve", "ILX:0794398"),
    ("posterior B gastric branch of vagus nerve", "ILX:0794399"),
    ("posterior C gastric branch of vagus nerve", "ILX:0794400"),
    # ===
    ("bronchial branch of right vagus nerve", "FMA:53613", "ILX:0791364"),
    ("superior cervical cardiac branch of right vagus nerve", "FMA:53598", "ILX:0786396"),
    ("right recurrent laryngeal nerve", "UBERON:0011767", "ILX:0728322"),
    ("extra laryngeal branch of right recurrent laryngeal nerve to larynx", "ILX:0795068"),
    ("branch of right recurrent laryngeal nerve to muscle of larynx", "ILX:0795084"),
    ("esophageal branch of right recurrent laryngeal nerve", "FMA:53608", "ILX:0785794"),
    ("tracheal branch of right recurrent laryngeal nerve", "FMA:53610", "ILX:0792888"),
    ("right inferior laryngeal nerve", "UBERON:0011767"),  # same code as for left recurrent laryngeal nerve
    ("external branch of right inferior laryngeal nerve", "None"),  # NA removed
    ("anterior branch of right recurrent laryngeal nerve", "ILX:0795070"),
    ("posterior branch of right recurrent laryngeal nerve", "ILX:0795072"),
    ("communicating branch of right recurrent laryngeal nerve with superior cervical ganglion", "None"),  # NA rejected
    ("communicating branch of right recurrent laryngeal nerve with right internal laryngeal nerve", "FMA:53528", "ILX:0790555"),
    ("inferior cervical cardiac branch of right recurrent laryngeal nerve", "ILX:0794235"),
    ("inferior cervical cardiac branch of right vagus nerve", "FMA:6713", "ILX:0789456"),
    ("thoracic cardiac branch of right vagus nerve", "FMA:53604", "ILX:0790057"),
    ("cardiac branch of right vagus to deep cardiac plexus", "FMA:6711", "ILX:0791784"),
    ("right branch of right vagus nerve to esophageal nerve plexus", "ILX:0794299"),

    # posterior vagus, a continuation of right vagus
    ("posterior esophageal vagus trunk", "ILX:0794858"),
    ("celiac branch of posterior vagal trunk", "FMA:6667", 'ILX:0789580'),
    ("greater posterior gastric nerve", "FMA:6689", "ILX:0788809"),
    ("pyloric branch of greater posterior gastric nerve", "FMA:6677"),

    # left vagus branches
    ("left vagus nerve", "FMA:6220", "ILX:0785628"),
    ("left vagus X nerve trunk", "UBERON:0035020"),
    ("left cervical trunk", "ILX:0794142"),
    ("left thoracic trunk", "ILX:0787543"),
    ("left meningeal branch of left vagus nerve", "FMA:53542", "ILX:0736691"),
    ("left branch between vagus nerve and glossopharyngeal nerve", "FMA:53560", "ILX:0790685"),
    ("left auricular branch of left vagus nerve", "FMA:53535", "ILX:0789344"),
    ("communicating branch of auricular branch of left vagus nerve with left facial nerve", "FMA:53588", "ILX:0787673"),
    ("communicating branch of auricular branch of left vagus nerve with left posterior auricular nerve", "FMA:53590", "ILX:0786219"),
    ("left pharyngeal branch of left vagus nerve to pharyngeal nerve plexus", "FMA:53636", "ILX:0789210"),
    ("lingual branch of left vagus nerve", "FMA:53634", "ILX:0791548"),
    ("left pharyngeal branch of left vagus nerve to superior cervical ganglion", "ILX:0795067"),
    ("branch of left vagus nerve to carotid body", "FMA:53607", "ILX:0791085"),
    ("left branch between vagus nerve and superior cervical ganglion", "ILX:0794061"),
    ("left superior laryngeal nerve", "FMA:53536", "ILX:0788780"),
    ("left internal laryngeal nerve", "FMA:53540", "ILX:0791167"),
    ("left external laryngeal nerve", "FMA:53538", "ILX:0789760"),
    ("superior branch of left internal laryngeal nerve", "FMA:53576", "ILX:0785786"),
    ("middle branch of left internal laryngeal nerve", "None"),  # NA pending
    ("inferior branch of left internal laryngeal nerve", "FMA:53582", "ILX:0785467"),
    ("communicating branch of left internal laryngeal nerve with left recurrent laryngeal nerve", "FMA:53572", "ILX:0789900"),
    ("communicating branch of left external laryngeal nerve with left superior cardiac nerve", "FMA:53562", "ILX:0787107"),
    ("left A cervical cardiopulmonary branch of vagus nerve", "ILX:0794160"),
    ("left B cervical cardiopulmonary branch of vagus nerve", "ILX:0794161"),
    ("pulmonary branch of left vagus nerve", "FMA:6679", "ILX:0792971"),
    ("left pulmonary branch A of the vagus nerve", "ILX:0795079"),
    ("left pulmonary branch B of the vagus nerve", "ILX:0795080"),
    ("left pulmonary branch C of the vagus nerve", "ILX:0795081"),
    ("left pulmonary branch D of the vagus nerve", "ILX:0795082"),
    ("left pulmonary branch E of the vagus nerve", "ILX:0795083"),
    # === temporary branches used in the japanese dataset
    ("left A medial pulmonary branch of vagus nerve", "ILX:0794266"),
    ("left B medial pulmonary branch of vagus nerve", "ILX:0794267"),
    # ===
    ("bronchial branch of left vagus nerve", "FMA:53558", "ILX:0791685"),
    ("superior cervical cardiac branch of left vagus nerve", "FMA:53599", "ILX:0791998"),
    ("left recurrent laryngeal nerve", "UBERON:0011766", "ILX:0724431"),
    ("extra laryngeal branch of left recurrent laryngeal nerve to larynx", "ILX:0795069"),
    ("branch of left recurrent laryngeal nerve to muscle of larynx", "ILX:0795085"),
    ("esophageal branch of left recurrent laryngeal nerve", "FMA:53609", "ILX:0791513"),
    ("tracheal branch of left recurrent laryngeal nerve", "FMA:53611", "ILX:0785091"),
    ("left inferior laryngeal nerve", "UBERON:0011766"),
    ("external branch of left inferior laryngeal nerve", "None"),  # NA removed
    ("anterior branch of left recurrent laryngeal nerve", "ILX:0795071"),
    ("posterior branch of left recurrent laryngeal nerve", "ILX:0795073"),
    ("communicating branch of left recurrent laryngeal nerve with superior cervical ganglion", "None"),  # NA rejected
    ("communicating branch of left recurrent laryngeal nerve with left internal laryngeal nerve", "FMA:53529", "ILX:0790440"),
    ("inferior cervical cardiac branch of left recurrent laryngeal nerve", "ILX:0794244"),
    ("inferior cervical cardiac branch of left vagus nerve", "FMA:6714", "ILX:0786047"),
    ("thoracic cardiac branch of left vagus nerve", "FMA:53605", "ILX:0791489"),
    ("cardiac branch of left vagus to deep cardiac plexus", "ILX:0795087"),
    ("left branch of left vagus nerve to esophageal nerve plexus", "ILX:0794310"),

    # anterior vagus, a continuation of left vagus
    ("anterior esophageal vagus trunk", "ILX:0794854"),
    ("hepatic branch of anterior vagal trunk", "FMA:6666", "ILX:0784595"),
    ("greater anterior gastric nerve", "FMA:6684", "ILX:0793831"),
    ("branch of greater anterior gastric nerve to coeliac nerve plexus", "FMA:53675"),

]


def get_vagus_branch_term(name):
    """
    Find term by matching name to any identifier held for a term.
    return: ( preferred name, preferred id )
    """
    for term in vagus_branch_terms:
        if name in term:
            return term[0], term[1]
    return name, ""


def get_vagus_marker_term(name: str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    return: ( preferred name, preferred id )
    """
    for term in vagus_marker_terms:
        if name in term:
            return term[0], term[1]
    raise NameError("Vagus annotation term '" + name + "' not found.")


def marker_name_in_terms(name: str):
    """
    Check if term exists in approved marker terms
    """
    for term in vagus_marker_terms:
        if name in term:
            return True
    return False



