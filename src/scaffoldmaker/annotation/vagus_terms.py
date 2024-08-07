"""
Common resource for vagus annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
vagus_terms = [
    # anatomical landmarks
    ("level of superior border of jugular foramen on the vagus nerve", "ILX:0794617"),
    ("level of inferior border of jugular foramen on the vagus nerve", "ILX:0794620"),
    ("level of inferior border of cranium on the vagus nerve", "ILX:0794623"),
    ("level of C1 transverse process on the vagus nerve", "ILX:0794626"),
    ("level of greater horn of hyoid on the vagus nerve", "ILX:0794629"),
    ("level of laryngeal prominence on the vagus nerve", "ILX:0794632"),
    ("level of angle of the mandible on the vagus nerve", "ILX:0794635"),
    ("level of carotid bifurcation on the vagus nerve", "ILX:0794638"),
    ("level of superior border of the clavicle on the vagus nerve", "ILX:0794641"),
    ("level of jugular notch on the vagus nerve", "ILX:0794644"),
    ("level of carina on the vagus nerve", ""),
    ("level of sternal angle on the vagus nerve", "ILX:0794647"),
    ("1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794650"),  # !!!
    ("level of esophageal hiatus on the vagus nerve", "ILX:0794653"),
    ("level of aortic hiatus on the vagus nerve", "ILX:0794656"),

    ("right level of superior border of jugular foramen on the vagus nerve", "ILX:0794618"),
    ("right level of inferior border of jugular foramen on the vagus nerve", "ILX:0794621"),
    ("right level of inferior border of cranium on the vagus nerve", "ILX:0794624"),
    ("right level of C1 transverse process on the vagus nerve", "ILX:0794627"),
    ("right level of greater horn of hyoid on the vagus nerve", "ILX:0794630"),
    ("right level of laryngeal prominence on the vagus nerve", "ILX:0794633"),
    ("right level of angle of the mandible on the vagus nerve", "ILX:0794636"),
    ("right level of carotid bifurcation on the vagus nerve", "ILX:0794639"),
    ("right level of superior border of the clavicle on the vagus nerve", "ILX:0794642"),
    ("right level of jugular notch on the vagus nerve", "ILX:0794645"),
    ("right level of carina on the vagus nerve", ""),
    ("right level of sternal angle on the vagus nerve", "ILX:0794648"),
    ("right 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794651"),  # !!!
    ("right level of esophageal hiatus on the vagus nerve", "ILX:0794654"),
    ("right level of aortic hiatus on the vagus nerve", "ILX:0794657"),

    ("left level of superior border of jugular foramen on the vagus nerve", "ILX:0794619"),
    ("left level of inferior border of jugular foramen on the vagus nerve", "ILX:0794622"),
    ("left level of inferior border of cranium on the vagus nerve", "ILX:0794625"),
    ("left level of C1 transverse process on the vagus nerve", "ILX:0794628"),
    ("left level of greater horn of hyoid on the vagus nerve", "ILX:0794631"),
    ("left level of laryngeal prominence on the vagus nerve","ILX:0794634"),
    ("left level of angle of the mandible on the vagus nerve", "ILX:0794637"),
    ("left level of carotid bifurcation on the vagus nerve", "ILX:0794640"),
    ("left level of superior border of the clavicle on the vagus nerve", "ILX:0794643"),
    ("left level of jugular notch on the vagus nerve", "ILX:0794646"),
    ("left level of carina on the vagus nerve", ""),
    ("left level of sternal angle on the vagus nerve", "ILX:0794649"),
    ("left 1 cm superior to start of esophageal plexus on the vagus nerve", "ILX:0794652"), # !!!
    ("left level of esophageal hiatus on the vagus nerve", "ILX:0794655"),
    ("left level of aortic hiatus on the vagus nerve", "ILX:0794658"),

    # branch names, based on lit.review
    # branch names with no indication of the side
    ("Meningeal branch of vagus nerve", "FMA:6231"),
    ("Communicating branch of vagus nerve with glossopharyngeal nerve", "FMA:6233"),
    ("Auricular branch of vagus nerve", "FMA:6232"),
    ("Pharyngeal branch of vagus nerve", "UBERON:0000929"),
    ("Lingual branch of vagus nerve", "FMA:6235"),
    ("Branch of vagus nerve to carotid body", "FMA:6237"),
    ("Communicating branch of superior cervical ganglion with vagus nerve", "FMA:6901"),
    ("Superior laryngeal nerve", "UBERON:0011326"),
    ("Internal laryngeal nerve", "FMA:6240"),
    ("External branch of superior laryngeal nerve", "FMA:6243"),
    ("Communicating branch of internal laryngeal nerve with recurrent laryngeal nerve", "FMA:53544"),
    ("Communicating branch of external laryngeal nerve with superior cardiac nerve", "FMA:6708"),
    ("Pulmonary branch of vagus nerve", "FMA:65515"),
    ("Bronchial branch of vagus nerve", "FMA:6247"),
    ("Superior cervical cardiac branch of vagus nerve", "FMA:6244"),
    ("Recurrent laryngeal nerve", "UBERON:0003716"),
    ("Extra-laryngeal branch of recurrent laryngeal nerve to larynx", "FMA:6710"),
    ("Oesophageal branch of recurrent laryngeal nerve", "FMA:6248"),
    ("Tracheal branch of recurrent laryngeal nerve", "FMA:6249"),
    ("Inferior laryngeal nerve", "ILX:0776141"), # UBERON:0003716 (synonym for recurrent laryngeal nerve)
    ("Communicating branch of recurrent laryngeal nerve with internal laryngeal nerve", "FMA:53526"),
    ("Inferior cervical cardiac branch of vagus nerve", "FMA:75530"),
    ("Thoracic cardiac branch of vagus nerve", "FMA:53601"),

    # right vagus
    ("Right vagus nerve", "FMA:6219"),
    ("right vagus X nerve trunk", "UBERON:0035021"),
    ("Meningeal branch of right vagus nerve", "FMA:53541"),
    ("Recurrent meningeal branch of right vagus nerve", ""),
    ("Communicating branch of right vagus nerve with right glossopharyngeal nerve", "FMA:52559"),
    ("Auricular branch of right vagus nerve", "FMA:53534"),
    ("Communicating branch of auricular branch of right vagus nerve with right facial nerve", ""),
    ("Communicating branch of auricular branch of right vagus nerve with right posterior auricular nerve", ""),
    ("Pharyngeal branch of right vagus nerve to pharyngeal nerve plexus", "FMA:53635"),
    ("Lingual branch of right vagus nerve", "FMA:53633"),
    ("Communicating branch of pharyngeal branch of right vagus nerve with superior cervical ganglion", ""),
    ("Communicating branch of pharyngeal branch of right vagus nerve", ""),
    ("Branch of right vagus nerve to carotid body", "FMA:53606"),
    ("Communicating branch of superior cervical ganglion with right vagus nerve", ""),
    ("Right superior laryngeal nerve", "FMA:53530"),
    ("Right internal laryngeal nerve", "FMA:53539"),
    ("Right external laryngeal nerve", "FMA:53537"),
    ("Superior ramus of right internal laryngeal nerve", ""),
    ("Middle ramus of right internal laryngeal nerve", ""),
    ("Inferior ramus of right internal laryngeal nerve", ""),
    ("Communicating branch of right internal laryngeal nerve with right recurrent laryngeal nerve", "FMA:53571"),
    ("Communicating branch of right external laryngeal nerve with right superior cardiac nerve", "FMA:53561"),
    ("Cardiopulmonary branch A of right vagus nerve", ""),
    ("Cardiopulmonary branch B of right vagus nerve", ""),
    ("Pulmonary branch of right vagus nerve", "FMA:6671"),
    ("Pulmonary branch A of right vagus nerve", ""),
    ("Pulmonary branch B of right vagus nerve", ""),
    ("Pulmonary branch C of right vagus nerve", ""),
    ("Pulmonary branch D of right vagus nerve", ""),
    ("Pulmonary branch E of right vagus nerve", ""),
    # === temporary branches used in the japanese dataset
    ("right A lateral pulmonary branch of vagus nerve", "ILX:0794278"),
    ("right B lateral pulmonary branch of vagus nerve", "ILX:0794279"),
    ("posterior A gastric branch of vagus nerve", "ILX:0794398"),
    ("posterior B gastric branch of vagus nerve", "ILX:0794399"),
    ("posterior C gastric branch of vagus nerve", "ILX:0794400"),
    # ===
    ("Bronchial branch of right vagus nerve", "FMA:53613"),
    ("Superior cervical cardiac branch of right vagus nerve", "FMA:53598"),
    ("Right recurrent laryngeal nerve", "UBERON:0011767"),
    ("Extra-laryngeal branch of right recurrent laryngeal nerve to larynx", "FMA:6250"), # code with no side indication
    ("Branch of right recurrent laryngeal nerve to muscle of larynx", ""),
    ("Oesophageal branch of right recurrent laryngeal nerve", "FMA:53608"),
    ("Tracheal branch of right recurrent laryngeal nerve", "FMA:53610"),
    # No independent code - same code as for left recurrent laryngeal nerve, for now used a different one
    ("Right inferior laryngeal nerve", "ILX:0776141"),
    ("External branch of right inferior laryngeal nerve", ""),
    ("Anterior branch of right inferior laryngeal nerve", ""),
    ("Posterior branch of right inferior laryngeal nerve", ""),
    ("Communicating branch of right recurrent laryngeal nerve with superior cervical ganglion", ""),
    ("Communicating branch of right recurrent laryngeal nerve with right internal laryngeal nerve", "FMA:53528"),
    ("Inferior cervical cardiac branch of right recurrent laryngeal nerve", ""),
    ("Inferior cervical cardiac branch of right vagus nerve", "FMA:6713"),
    ("Thoracic cardiac branch of right vagus nerve", "FMA:53604"),
    ("Right cardiac branch to deep cardiac nerve plexus", "ILX:0794165"), # not sure about id
    ("Branch of right vagus nerve to oesophageal nerve plexus", ""),

    # left vagus
    ("Left vagus nerve", "FMA:6220"),
    ("left vagus X nerve trunk", "UBERON:0035020"),
    ("Meningeal branch of left vagus nerve", "FMA:53542"),
    ("Recurrent meningeal branch of left vagus nerve", ""),
    ("Communicating branch of left vagus nerve with left glossopharyngeal nerve", "FMA:53560"),
    ("Auricular branch of left vagus nerve", "FMA:53535"),
    ("Communicating branch of auricular branch of left vagus nerve with left facial nerve", ""),
    ("Communicating branch of auricular branch of left vagus nerve with left posterior auricular nerve", ""),
    ("Pharyngeal branch of left vagus nerve to pharyngeal nerve plexus", "FMA:53636"),
    ("Lingual branch of left vagus nerve", "FMA:53634"),
    ("Communicating branch of pharyngeal branch of left vagus nerve with superior cervical ganglion", ""),
    ("Communicating branch of pharyngeal branch of left vagus nerve", ""),
    ("Branch of left vagus nerve to carotid body", "FMA:53607"),
    ("Communicating branch of superior cervical ganglion with left vagus nerve", ""),
    ("Left superior laryngeal nerve", "FMA:53536"),
    ("Left internal laryngeal nerve", "FMA:53540"),
    ("Left external laryngeal nerve", "FMA:53538"),
    ("Superior ramus of left internal laryngeal nerve", ""),
    ("Middle ramus of left internal laryngeal nerve", ""),
    ("Inferior ramus of left internal laryngeal nerve", ""),
    ("Communicating branch of left internal laryngeal nerve with left recurrent laryngeal nerve", "FMA:53572"),
    ("Communicating branch of left external laryngeal nerve with left superior cardiac nerve", "FMA:53562"),
    ("Cardiopulmonary branch A of left vagus nerve", ""),
    ("Cardiopulmonary branch B of left vagus nerve", ""),
    ("Pulmonary branch of left vagus nerve", "FMA:6679"),
    ("Pulmonary branch A of left vagus nerve", ""),
    ("Pulmonary branch B of left vagus nerve", ""),
    ("Pulmonary branch C of left vagus nerve", ""),
    ("Pulmonary branch D of left vagus nerve", ""),
    # === temporary branches used in the japanese dataset
    ("left A medial pulmonary branch of vagus nerve", "ILX:0794266"),
    ("left B medial pulmonary branch of vagus nerve", "ILX:0794267"),
    # ===
    ("Bronchial branch of left vagus nerve", "FMA:53558"),
    ("Superior cervical cardiac branch of left vagus nerve", "FMA:53599"),
    ("Left recurrent laryngeal nerve", "UBERON:0011766"),
    ("Extra-laryngeal branch of left recurrent laryngeal nerve to larynx", ""),
    ("Branch of left recurrent laryngeal nerve to muscle of larynx", "FMA:6250"), # code with no side indication
    ("Oesophageal branch of left recurrent laryngeal nerve", "FMA:53609"),
    ("Tracheal branch of left recurrent laryngeal nerve", "FMA:53611"),
    # No independent code - same code as for left recurrent laryngeal nerve, for now used a different one
    ("Left inferior laryngeal nerve", "ILX:0776141"),
    ("External branch of left inferior laryngeal nerve", ""),
    ("Anterior branch of left inferior laryngeal nerve", ""),
    ("Posterior branch of left inferior laryngeal nerve", ""),
    ("Communicating branch of left recurrent laryngeal nerve with superior cervical ganglion", ""),
    ("Communicating branch of left recurrent laryngeal nerve with left internal laryngeal nerve", "FMA:53529"),
    ("Inferior cervical cardiac branch of left recurrent laryngeal nerve", ""),
    ("Inferior cervical cardiac branch of left vagus nerve", "FMA:6714"),
    ("Thoracic cardiac branch of left vagus nerve", "FMA:53605"),
    ("Left cardiac branch to deep cardiac nerve plexus", "ILX:0794171"), # not sure about id
    ("Branch of left vagus nerve to oesophageal nerve plexus", ""),

    # anterior vagus, a continuation of left vagus
    ("Anterior esophageal vagus trunk", "ILX:0794854"),
    ("Hepatic branch of anterior vagal trunk", "FMA:6666"),
    ("Greater anterior gastric nerve", "FMA:6684"),
    ("Branch of greater anterior gastric nerve to coeliac nerve plexus", "FMA:53675"),

    # posterior vagus, a continuation of right vagus
    ("Posterior esophageal vagus trunk", "ILX:0794858"),
    ("Coeliac branch of posterior vagal trunk", "FMA:6667"),
    ("Greater posterior gastric nerve", "FMA:6689"),
    ("Pyloric branch of greater posterior gastric nerve", "FMA:6677"),
]


def get_vagus_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in vagus_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Vagus annotation term '" + name + "' not found.")