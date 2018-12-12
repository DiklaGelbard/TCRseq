from collections import Counter, defaultdict

import six
import pandas as pd

class Basic_Cell(object):
    def __init__(self, cell_name, barcode_df, species="Hsap"):
        self.name = cell_name
        self.species = species
        self.barcode_df = barcode_df


class Cell(object):
    """Class to describe T cells containing A and B loci"""
    def __init__(self, cell_name, recombinants, is_empty=False, species="Mmus",
                 receptor=None, loci=None):

        self.name = cell_name
        self.recombinants = self._process_recombinants(recombinants, receptor,
                                                       loci)
        self.is_empty = self._check_is_empty()
        self.species = species

    def _process_recombinants(self, recombinants, receptor, loci):
        recombinant_dict = defaultdict(dict)
        if recombinants is not None:
            for r_name, r in six.iteritems(recombinants):
                r_name = r_name.split("_")
                receptor = r_name[0]
                locus = r_name[1]
                recombinant_dict[receptor][locus] = r

        # normalise this to put None in cases where no receptors found
        for l in loci:
            if l not in recombinant_dict[receptor]:
                recombinant_dict[receptor][l] = None
        return dict(recombinant_dict)

    def _check_is_empty(self):
        if (self.recombinants is None or len(self.recombinants) == 0):
            return True
        else:
            return False

    def __str__(self):
        return (self.name)

    def choose_recombinants(self):
        ret_dict = defaultdict(dict)
        for receptor, locus_dict in six.iteritems(self.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None and len(recombinants) > 0:
                    V_ranks = Counter()
                    D_ranks = Counter()
                    J_ranks = Counter()
                    cdr3_ranks = Counter()
                    cdr3_unique_reads = Counter()
                    cdr3_percent_identity = Counter()
                    cdr3_translation_ranks = Counter()
                    cdr3_translation_unique_reads = Counter()
                    cdr3_translation_percent_identity = Counter()
                    V_eValue = Counter()
                    D_eValue = Counter()
                    J_eValue = Counter()
                    for rec in recombinants:
                        ## query name q contig name
                        freq=rec.contig_name.strip().split("-")[1]
                        for hit in rec.hit_table:
                            seq = hit[2]
                            if hit[0] == 'V':
                                V_ranks.update({str(seq): int(freq)})
                                V_eValue.update({str(seq): float(hit[12])*int(freq)})
                            elif hit[0] == 'D':
                                D_ranks.update({str(seq): int(freq)})
                                D_eValue.update({str(seq): float(hit[12])*int(freq)})
                            elif hit[0] == 'J':
                                J_ranks.update({str(seq):  int(freq)})
                                J_eValue.update({str(seq): float(hit[12])*int(freq)})
                        if rec.cdr3 is not None and len(rec.cdr3) != 0:
                            cdr3_ranks.update({rec.cdr3[1]: int(freq)})
                            cdr3_unique_reads.update({rec.cdr3[1]: 1})
                            cdr3_translation_ranks.update({rec.cdr3[2]: int(freq)})
                            cdr3_translation_unique_reads.update({rec.cdr3[2]: 1})
                            for align in rec.alignment_summary:
                                if 'CDR3' in align:
                                    a = align.split('\t')
                                    cdr3_percent_identity.update({rec.cdr3[1]: float(a[7])*int(freq)})
                                    cdr3_translation_percent_identity.update({rec.cdr3[2]: float(a[7])*int(freq)})
                    v_avg = dict()
                    for v in V_ranks.keys():
                        v_avg[v] = float(V_eValue[v]/V_ranks[v])
                    d_avg = dict()
                    for d in D_ranks.keys():
                        d_avg[d] = float(D_eValue[d]/D_ranks[d])
                    j_avg = dict()
                    for j in J_ranks.keys():
                        j_avg[j] = float(J_eValue[j]/J_ranks[j])
                    cdr3_avg_identity = dict()
                    cdr3_translation_avg_identity = dict()
                    for cdr3 in cdr3_ranks.keys():
                        cdr3_avg_identity[cdr3] = float(cdr3_percent_identity[cdr3]/cdr3_ranks[cdr3])
                    for cdr3 in cdr3_translation_ranks.keys():
                        cdr3_translation_avg_identity[cdr3] = float(cdr3_translation_percent_identity[cdr3]/cdr3_translation_ranks[cdr3])
                    ret_dict[receptor][locus] = pd.DataFrame([{"cell_name": self.name ,
                                                               "V_first": "NA" if len(V_ranks) == 0 else V_ranks.most_common(2)[0][0],
                                                               "V_first_counts": 0 if len(V_ranks) == 0  else V_ranks.most_common(2)[0][1],
                                                               "V_first_avg_e_value": 0 if len(V_ranks) == 0  else v_avg[V_ranks.most_common(2)[0][0]],
                                                               "V_second": "NA" if len(V_ranks) < 2 else V_ranks.most_common(2)[1][0],
                                                               "V_second_counts": 0 if len(V_ranks) < 2 else  V_ranks.most_common(2)[1][1],
                                                               "V_second_avg_e_value": 0 if len(V_ranks) < 2  else v_avg[V_ranks.most_common(2)[1][0]],
                                                               "D_first": "NA" if len(D_ranks) == 0 else D_ranks.most_common(2)[0][0],
                                                               "D_first_counts": 0 if len(D_ranks) == 0 else D_ranks.most_common(2)[0][1],
                                                               "D_first_avg_e_value": 0 if len(D_ranks) == 0  else d_avg[D_ranks.most_common(2)[0][0]],
                                                               "D_second": "NA" if len(D_ranks) < 2 else  D_ranks.most_common(2)[1][0],
                                                               "D_second_counts": 0 if len(D_ranks) < 2 else D_ranks.most_common(2)[1][1],
                                                               "D_second_avg_e_value": 0 if len(D_ranks) < 2  else d_avg[D_ranks.most_common(2)[1][0]],
                                                               "J_first": "NA"  if len(J_ranks) == 0 else J_ranks.most_common(2)[0][0],
                                                               "J_first_counts": 0 if len(J_ranks) == 0 else J_ranks.most_common(2)[0][1],
                                                               "J_first_avg_e_value": 0 if len(J_ranks) == 0  else j_avg[J_ranks.most_common(2)[0][0]],
                                                               "J_second": "NA"  if len(J_ranks) < 2 else J_ranks.most_common(2)[1][0],
                                                               "J_second_counts": 0  if len(J_ranks) < 2 else J_ranks.most_common(2)[1][1],
                                                               "J_second_avg_e_value": 0 if len(J_ranks) < 2  else j_avg[J_ranks.most_common(2)[1][0]],
                                                               "CDR3_first": "NA" if len(cdr3_ranks) == 0 else cdr3_ranks.most_common(2)[0][0],
                                                               "CDR3_first_counts": 0 if len(cdr3_ranks) == 0 else cdr3_ranks.most_common(2)[0][1],
                                                               "CDR3_first_identity": 0 if len(cdr3_ranks) == 0 else cdr3_avg_identity[cdr3_ranks.most_common(2)[0][0]],
                                                               "CDR3_first_unique_reads": 0 if len(cdr3_ranks) == 0 else cdr3_unique_reads[cdr3_ranks.most_common(2)[0][0]],
                                                               "CDR3_translation_first": "NA" if len(cdr3_translation_ranks) == 0 else
                                                               cdr3_translation_ranks.most_common(2)[0][0],
                                                               "CDR3_translation_first_counts": 0 if len(cdr3_translation_ranks) == 0 else
                                                               cdr3_translation_ranks.most_common(2)[0][1],
                                                               "CDR3_translation_first_identity": 0 if len(cdr3_translation_ranks) == 0 else
                                                               cdr3_translation_avg_identity[cdr3_translation_ranks.most_common(2)[0][0]],
                                                               "CDR3_translation_first_unique_reads": 0 if len(
                                                                   cdr3_translation_ranks) == 0 else cdr3_translation_unique_reads[
                                                                   cdr3_translation_ranks.most_common(2)[0][0]],
                                                               "CDR3_second":"NA" if len(cdr3_ranks) < 2 else cdr3_ranks.most_common(2)[1][0],
                                                               "CDR3_second_counts": 0 if len(cdr3_ranks) < 2  else cdr3_ranks.most_common(2)[1][1],
                                                               "CDR3_second_identity": 0 if len(cdr3_ranks) < 2 else cdr3_avg_identity[cdr3_ranks.most_common(2)[1][0]],
                                                               "CDR3_second_unique_reads": 0 if len(cdr3_ranks) < 2 else cdr3_unique_reads[cdr3_ranks.most_common(2)[1][0]],
                                                               "CDR3_translation_second": "NA" if len(cdr3_translation_ranks) < 2 else
                                                               cdr3_translation_ranks.most_common(2)[1][0],
                                                               "CDR3_translation_second_counts": 0 if len(cdr3_translation_ranks) < 2  else
                                                               cdr3_translation_ranks.most_common(2)[1][1],
                                                               "CDR3_translation_second_identity": 0 if len(cdr3_translation_ranks) < 2 else
                                                               cdr3_translation_avg_identity[cdr3_translation_ranks.most_common(2)[1][0]],
                                                               "CDR3_translation_second_unique_reads": 0 if len(
                                                                   cdr3_translation_ranks) < 2 else cdr3_translation_unique_reads[
                                                                   cdr3_translation_ranks.most_common(2)[1][0]]}],
                                                             columns=["cell_name","V_first","V_first_counts","V_first_avg_e_value",
                                                                      "V_second","V_second_counts","V_second_avg_e_value",
                                                                      "D_first","D_first_counts","D_first_avg_e_value",
                                                                      "D_second","D_second_counts","D_second_avg_e_value",
                                                                      "J_first", "J_first_counts","J_first_avg_e_value",
                                                                      "J_second", "J_second_counts","J_second_avg_e_value",
                                                                      "CDR3_first", "CDR3_first_counts","CDR3_first_identity","CDR3_first_unique_reads",
                                                                      "CDR3_translation_first","CDR3_translation_first_counts","CDR3_translation_first_identity","CDR3_translation_first_unique_reads",
                                                                      "CDR3_second", "CDR3_second_counts","CDR3_second_identity","CDR3_second_unique_reads",
                                                                      "CDR3_translation_second","CDR3_translation_second_counts","CDR3_translation_second_identity",
                                                                      "CDR3_translation_second_unique_reads"])
        print(ret_dict)
        return ret_dict


class Recombinant(object):
    """Class to describe a recombined TCR locus as determined from the single-cell pipeline"""

    def __init__(self, contig_name, locus, identifier, all_poss_identifiers,
                 productive, stop_codon, in_frame, TPM,
                 dna_seq, hit_table, summary, junction_details, best_VJ_names,
                 alignment_summary, fasta_seq,
                 imgt_reconstructed_seq, has_D,cdr3):
        self.contig_name = contig_name
        self.locus = locus
        self.identifier = identifier
        self.all_poss_identifiers = all_poss_identifiers
        self.productive = productive
        self.TPM = TPM
        self.dna_seq = dna_seq
        self.cdr3 = cdr3
        self.hit_table = hit_table
        self.summary = summary
        self.junction_details = junction_details
        self.best_VJ_names = best_VJ_names
        self.alignment_summary = alignment_summary
        self.in_frame = in_frame
        self.stop_codon = stop_codon
        self.fasta_seq = fasta_seq
        self.imgt_reconstructed_seq = imgt_reconstructed_seq
        self.has_D_segment = has_D

    def __str__(self):
        return (
        "{} {} {} {}".format(self.identifier, self.productive, self.TPM))


    def get_summary(self):
        summary_string = "##{contig_name}##\n".format(
            contig_name=self.contig_name)
        if not self.has_D_segment:
            V_segment = self.summary[0]
            J_segment = self.summary[1]
            segments_string = "V segment:\t{V_segment}\n" \
                              "J segment:\t{J_segment}\n".format(
                V_segment=V_segment,
                J_segment=J_segment)
        else:
            V_segment = self.summary[0]
            D_segment = self.summary[1]
            J_segment = self.summary[2]
            segments_string = "V segment:\t{V_segment}\nD segment:\t{D_segment}\n" \
                              "J segment:\t{J_segment}\n".format(
                V_segment=V_segment, D_segment=D_segment, J_segment=J_segment)
        summary_string += segments_string
        summary_string += "dna_seq:\t{dna_seq}\nOriginal fasta seq:\t{fasta_seq}\nCDR3:\t{cdr3}\n".format(dna_seq=str(self.dna_seq).upper(),
                                                                                                          fasta_seq=self.fasta_seq,cdr3=self.cdr3)
        summary_string += "ID:\t{}\n".format(self.identifier)
        summary_string += "TPM:\t{TPM}\nProductive:\t{productive}\nStop codon:" \
                          "\t{stop_codon}\nIn frame:\t{in_frame}\n\n".format(
            TPM=self.TPM, productive=self.productive,
            stop_codon=self.stop_codon, in_frame=self.in_frame)

        summary_string += 'Segment\tquery_id\tsubject_id\t% identity\t' \
                          'alignment length\tmismatches\tgap opens\tgaps' \
                          '\tq start\tq end\ts start\ts end\te value\tbit score\n'
        for line in self.hit_table:
            summary_string = summary_string + "\t".join(line) + "\n"
        return (summary_string)


