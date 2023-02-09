import sbol3
from typing import List, Tuple
from sbol_utilities.component import dna_component_with_sequence
import tyto

def backbone(identity: str, sequence: str, dropout_location: List[int], fusion_site_length:int, linear:bool, **kwargs) -> Tuple[sbol3.Component, sbol3.Sequence]:
    """Creates a Backbone Component and its Sequence.

    :param identity: The identity of the Component. The identity of Sequence is also identity with the suffix '_seq'.
    :param sequence: The DNA sequence of the Component encoded in IUPAC.
    :param dropout_location: List of 2 integers that indicates the start and the end of the dropout sequence including overhangs. Note that the index of the first location is 1, as is typical practice in biology, rather than 0, as is typical practice in computer science.
    :param fusion_site_length: Integer of the lenght of the fusion sites (eg. BsaI fusion site lenght is 4, SapI fusion site lenght is 3)
    :param linear: Boolean than indicates if the backbone is linear, by default it is seted to Flase which means that it has a circular topology.
    :param kwargs: Keyword arguments of any other Component attribute.
    :return: A tuple of Component and Sequence.
    """
    if len(dropout_location) != 2:
        raise ValueError('The dropout_location only accepts 2 int values in a list.')
    backbone_component, backbone_seq = dna_component_with_sequence(identity, sequence, **kwargs)
    backbone_component.roles.append(sbol3.SO_DOUBLE_STRANDED)  
    dropout_location_comp = sbol3.Range(sequence=backbone_seq, start=dropout_location[0], end=dropout_location[1])
    insertion_site_location1 = sbol3.Range(sequence=backbone_seq, start=dropout_location[0], end=dropout_location[0]+fusion_site_length, order=1)
    insertion_site_location2 = sbol3.Range(sequence=backbone_seq, start=dropout_location[1]-fusion_site_length, end=dropout_location[1], order=3)
    dropout_sequence_feature = sbol3.SequenceFeature(locations=[dropout_location_comp], roles=[tyto.SO.deletion])
    insertion_sites_feature = sbol3.SequenceFeature(locations=[insertion_site_location1, insertion_site_location2], roles=[tyto.SO.insertion_site])
    if linear:
        backbone_component.types.append(sbol3.SO_LINEAR)
        backbone_component.roles.append(sbol3.SO_ENGINEERED_REGION)
        open_backbone_location1 = sbol3.Range(sequence=backbone_seq, start=1, end=dropout_location[0]+fusion_site_length-1, order=1)
        open_backbone_location2 = sbol3.Range(sequence=backbone_seq, start=dropout_location[1]-fusion_site_length, end=len(sequence), order=3)
        open_backbone_feature = sbol3.SequenceFeature(locations=[open_backbone_location1, open_backbone_location2])
    else: 
        backbone_component.types.append(sbol3.SO_CIRCULAR)
        backbone_component.roles.append(tyto.SO.plasmid_vector)
        open_backbone_location1 = sbol3.Range(sequence=backbone_seq, start=1, end=dropout_location[0]+fusion_site_length-1, order=2)
        open_backbone_location2 = sbol3.Range(sequence=backbone_seq, start=dropout_location[1]-fusion_site_length, end=len(sequence), order=1)
        open_backbone_feature = sbol3.SequenceFeature(locations=[open_backbone_location1, open_backbone_location2])
    backbone_component.features.append(dropout_sequence_feature)
    backbone_component.features.append(insertion_sites_feature)
    backbone_component.features.append(open_backbone_feature)
    backbone_dropout_meets = sbol3.Constraint(restriction='http://sbols.org/v3#meets', subject=dropout_sequence_feature, object=open_backbone_feature)
    backbone_component.constraints.append(backbone_dropout_meets)
    return backbone_component, backbone_seq

lvl1_pOdd_acceptor_seq = 'gctcgagtcccgtcaagtcagcgtaatgctctgccagtgttacaaccaattaaccaattctgattagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagcttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagcgagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgttttcccggggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaatcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctggagcaagacgtttcccgttgaatatggctcataacaccccttgtattactgtttatgtaagcagacagttttattgttcatgatgatatatttttatcttgtgcaatgtaacatcagagattttgagacacaacgtggctttgttgaataaatcgaacttttgctgagttgaaggatcagctcgagtgccacctgacgtctaagaaaccattattatcatgacattaacctataaaaataggcgtatcacgaggcagaatttcagataaaaaaaatccttagctttcgctaaggatgatttctggaattcgctcttcaatgggagtgagacccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacacatactagagaaagaggagaaatactagatggcttcctccgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaataacgctgatagtgctagtgtagatcgctactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataggtctcaCGCTgcatgaagagcctgcagtccggcaaaaaagggcaaggtgtcaccaccctgccctttttctttaaaaccgaaaagattacttcgcgttatgcaggcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccacaggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgaca'
podd_backbone, podd_backbone_seq = backbone('pOdd_bb', lvl1_pOdd_acceptor_seq, [1169,2259], 4, False, name='pOdd_bb')
