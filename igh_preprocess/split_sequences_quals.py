import sys

def main(argv):
    num_splits    = int(sys.argv[1])
    in_seq_quals = sys.argv[2]
    uids_split_on = sys.argv[3]
    
    outfiles = [open('%s.%s' % (in_seq_quals,n),'w') for n in range(num_splits)]
    
    uids = get_split_uids(uids_split_on)
    
    current_outfile = 0
    current_uid = 0 # index for which uid is expected next
    
    #Once a uid is found that is split after, continue writing lines to the output 
    #file until another uid is encountered. This is to handle the case where a uid has 
    #multiple quality scores
    write_continuation = 0
        
    with open(in_seq_quals,'rU') as f:
        for line in f:
                  
            if '>' in line[0] and len(line)<15 and write_continuation == 0:
            # quality scores can contain the character '>' which necessitates 
            # the above conditions for ensuring that the line is an uid line
                uid = line.strip().split('>')[1]
                if uid == uids[current_uid]:
                    write_continuation = 1
               
                outfiles[current_outfile].write(line)
               
            elif '>' in line[0] and len(line)<15 and write_continuation ==1:
                write_continuation = 0
                outfiles[current_outfile].close()
                current_outfile += 1
                
                if current_uid == num_splits-2:
                    #current uid does not change because this is the last split
                    # Ex: for splits = 5 there will only be 4 uids split on
                    # -2 is for zero indexing (uid 3 is the 4th)
                    pass
                else:
                    current_uid += 1
                
                outfiles[current_outfile].write(line)
            
            else:
                # write quality score
                outfiles[current_outfile].write(line)
                
    
def get_split_uids(infile):
    uids=[]
    with open(infile) as f:
        for line in f:
            uids.append(line.strip())
    return uids
    
    
if __name__ == "__main__":
    main(sys.argv)

    
