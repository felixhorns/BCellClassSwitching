import sys
import os

# Input
infile = sys.argv[1]  #expects: {dir}/consensus_complete.txt
variable_barcode_len = int(sys.argv[2])  # expects: boolean, 1 for variable barcode length

# Output
output_dir = os.path.dirname(infile) + "/"
out_molecule_ids = open(output_dir+'sequences_molecule_ids.fasta', 'w')
out_abundances = open(output_dir+'sequences_abundances.FASTA', 'w')
out_reads_per_molecule = open(output_dir+'sequences_reads_per_molecule.TXT', 'w')
out_quals = open(output_dir+'sequences_quals.TXT', 'w')

#####

length = 0
for line in open(infile, 'rU'): length += 1

molecule_ids = {}
quals = {}


def trim_quals(qual_list, location):
    '''
    :param qual_list
    :param location: beginning, end, or both
    :return: trimmed qual list
    '''
    if location == 'beginning':
        qual_list = [x[4:] for x in qual_list]
    elif location == 'end':
        qual_list = [x[:-4] for x in qual_list]
    elif location == 'both':
        qual_list = [x[4:-4] for x in qual_list]
    return qual_list


def update_seq_quals(existing_seq, seq, existing_quals, new_quals):
    '''
    handles trimming of quality scores based on the common sequence
    between the sequence being added and the one present in the dict
    :returns the final trimmed sequence common to both seq inputs, and
            the trimmed lists of quality scores
    '''

    if existing_seq == seq:
        # seqs are the same, no trimming quals
        final_seq = seq

    elif len(seq) < len(existing_seq):
        # take seq, trim existing_quals
        final_seq = seq

        # determine which end(s) of existing_quals to trim
        if seq[-40:] != existing_seq[-40:]:
            existing_quals = trim_quals(existing_quals, 'end')
        if seq[:40] != existing_seq[:40]:
            existing_quals = trim_quals(existing_quals, 'beginning')

    elif len(existing_seq) < len(seq):
        # take existing seq, trim seq quals
        final_seq = existing_seq

        # determine which end(s) of new_quals to trim
        if seq[-40:] != existing_seq[-40:]:
            new_quals = trim_quals(new_quals, 'end')
        if seq[:40] != existing_seq[:40]:
            new_quals = trim_quals(new_quals, 'beginning')

    elif len(existing_seq) == len(seq):
        # instance of right and left smaller barcode on separate molecules
        if seq[:40] == existing_seq[4:44]:
            # trim end of seq, (or beginning of existing seq)
            final_seq = seq[:-4]
            new_quals = trim_quals(new_quals, 'end')
            existing_quals = trim_quals(existing_quals, 'beginning')

        elif existing_seq[:40] == seq[4:44]:
            final_seq = existing_seq[:-4]
            new_quals = trim_quals(new_quals, 'beginning')
            existing_quals = trim_quals(existing_quals, 'end')

    # check quals length agreement
    set_new = set([len(x) for x in new_quals])
    set_exist = set([len(x) for x in existing_quals])
    assert len(set_new) == 1, 'New quals multiple lengths: %s' % set_new
    assert len(set_exist) == 1, 'Existing quals multiple lengths: %s' % set_exist
    assert set_new == set_exist, 'diff len: %s\n%s\n%s\n%s\n%s\n%s' % (set_new, set_exist,seq,existing_seq, new_quals,existing_quals)

    return final_seq, new_quals, existing_quals


# Main loop ------------------------------------

with open(infile, 'rU') as infile:

    line_num = 0

    if variable_barcode_len:
        print 'using variable length barcode'
        # due to variable barcode length, 4 possibilities for each seq
        while line_num < length:
            molecule_id = infile.readline().strip().split('@')[1]
            seq = infile.readline().strip()
            plus_line = infile.readline()
            #na = infile.readline().strip()
            #qual = ''
            qual = infile.readline().strip()            

            # some garbage seqs have strings of C's or G's of 20+
            if len(set(seq[-20:])) == 1:
                line_num += 4
                continue
            if len(set(seq[:20])) == 1:
                line_num += 4
                continue

            seq_4minus_left = seq[4:]
            seq_4minus_right = seq[:-4]
            seq_4minus_each = seq[4:-4]

            # seq lookup dict structure
            # key = sequence, value = tuple of MID list, quals list
            # {'seq': ( [molecular_ids], [quals] )

            # check if 4 seq possibilities are in dictionary
            if seq in molecule_ids:
                molecule_ids[seq][0].append(molecule_id)
                molecule_ids[seq][1].append(qual)

            elif seq_4minus_left in molecule_ids:
                molecule_ids[seq_4minus_left][0].append(molecule_id)
                molecule_ids[seq_4minus_left][1].append(qual[4:])

            elif seq_4minus_right in molecule_ids:
                molecule_ids[seq_4minus_right][0].append(molecule_id)
                molecule_ids[seq_4minus_right][1].append(qual[:-4])

            elif seq_4minus_each in molecule_ids:
                molecule_ids[seq_4minus_each][0].append(molecule_id)
                molecule_ids[seq_4minus_each][1].append(qual[4:-4])

            else:
                # otherwise add all 4
                molecule_ids[seq] = ( [molecule_id], [qual] )  # tuple
                molecule_ids[seq_4minus_left] = ( [molecule_id], [qual[4:]] )
                molecule_ids[seq_4minus_right] = ( [molecule_id], [qual[:-4]] )
                molecule_ids[seq_4minus_each] = ( [molecule_id], [qual[4:-4]] )

            line_num += 4

        # collapse molecule_ids dictionary
        out_seq_mid = {}  # key: unique id, val: set of molecular ids
        mid_lookup = {}     # key: molecular id, val: uid to which it belongs

        uid = 0
        for seq, mid_qual in molecule_ids.items():
            # mid_qual is a tuple of two lists: mids and quals
            state = 0
            for mid in mid_qual[0]:
                if mid in mid_lookup:
                    id_hit = mid
                    state = 1
                    break
            if state == 1 and len(mid_qual[0]) != 1:
                # compare seq to seq currently in dict
                # trim to overlapping seq, update quals
                existing_seq = out_seq_mid[mid_lookup[id_hit]][2]
                existing_quals = out_seq_mid[mid_lookup[id_hit]][1]
                new_quals = mid_qual[1]

                agreed_seq, new_quals, existing_quals = update_seq_quals(existing_seq, seq, existing_quals, new_quals)

                # update seq
                out_seq_mid[mid_lookup[id_hit]][2] = agreed_seq

                # update existing quals
                out_seq_mid[mid_lookup[id_hit]][1] = existing_quals

                # add mids and updated quals to existing UID
                for mid, qual in zip(mid_qual[0], new_quals):
                    if mid not in out_seq_mid[mid_lookup[id_hit]][0]:
                        out_seq_mid[mid_lookup[id_hit]][0].append(mid)
                        out_seq_mid[mid_lookup[id_hit]][1].append(qual)

            elif state == 0:
                # first instance of these molecular id(s)
                out_seq_mid[uid] = [mid_qual[0], mid_qual[1], seq]
                for mid in mid_qual[0]:
                    mid_lookup[mid] = uid
                uid += 1

        # write out UIDs
        for seq_id, uid in enumerate(out_seq_mid):

            out_molecule_ids.write(">"+str(seq_id)+"_"+",".join(out_seq_mid[uid][0])+"\n")
            out_molecule_ids.write(out_seq_mid[uid][2]+"\n")

            out_abundances.write(">" + str(seq_id) + "_" + str(len(out_seq_mid[uid][0])) + "\n")
            out_abundances.write(out_seq_mid[uid][2]+"\n")

            num_reads_per_molecule = [i.split("_")[-1] for i in out_seq_mid[uid][0]]
            out_reads_per_molecule.write(">" + str(seq_id) + "_" + str(len(out_seq_mid[uid][0])) + "\n")
            out_reads_per_molecule.write(",".join(num_reads_per_molecule) + "\n")

            out_quals.write(">"+str(seq_id)+"_"+str(len(out_seq_mid[uid][0]))+"\n")
            for qual in out_seq_mid[uid][1]: out_quals.write(qual+"\n")

    else:
        # barcodes are same length
        while line_num < length:
            molecule_id = infile.readline().strip().split('@')[1]
            seq = infile.readline().strip()
            plus_line = infile.readline()
            qual = infile.readline().strip()

            if seq in molecule_ids:
                molecule_ids[seq].append(molecule_id)
                quals[seq].append(qual)
            else:
                molecule_ids[seq] = [molecule_id]
                quals[seq] = [qual]

            line_num += 4

        for seq_id, seq in enumerate(molecule_ids.keys()):

            out_molecule_ids.write(">"+str(seq_id)+"_"+",".join(molecule_ids[seq])+"\n")
            out_molecule_ids.write(seq+"\n")

            out_abundances.write(">" + str(seq_id) + "_" + str(len(molecule_ids[seq])) + "\n")
            out_abundances.write(seq+"\n")

            num_reads_per_molecule = [i.split("_")[-1] for i in molecule_ids[seq]]
            out_reads_per_molecule.write(">" + str(seq_id) + "_" + str(len(molecule_ids[seq])) + "\n")
            out_reads_per_molecule.write(",".join(num_reads_per_molecule) + "\n")

            out_quals.write(">"+str(seq_id)+"_"+str(len(molecule_ids[seq]))+"\n")
            for qual in quals[seq]:
                out_quals.write(qual+"\n")


out_molecule_ids.close()
out_abundances.close()
out_reads_per_molecule.close()
out_quals.close()

