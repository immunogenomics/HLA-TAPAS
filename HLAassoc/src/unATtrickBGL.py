import os, sys, re
import gzip

normal_nts = ['A', 'C', 'G', 'T']

p_FristTwo = re.compile(r'^\S+\s+(\S+)\s+')



def unATtrickBGL(_bgl, _ref_bim, _out, _tar_fam=None):
    
    f_bgl = gzip.open(_bgl, 'rt') if _bgl.endswith('.gz') else open(_bgl, 'r')
    f_bim = open(_ref_bim, 'r')
    f_out = open(_out, 'w')
    
    flag_useFAM = bool(_tar_fam) and os.path.exists(_tar_fam)
    
    if flag_useFAM:
        
        l_FID = ['P', 'pedigree']
        l_IID = ['I', 'id']
        l_PID = ['f', 'fID']
        l_MID = ['m', 'mID']
        l_Sex = ['C', 'gender']
        
        with open(_tar_fam, 'r') as f_fam:
            
            for line in f_fam:
                
                [FID, IID, PID, MID, Sex, Phe] = re.split(r'\s+', line.rstrip('\n'))
                
                l_FID.append(FID); l_FID.append(FID);
                l_IID.append(IID); l_IID.append(IID);
                l_PID.append(PID); l_PID.append(PID);
                l_MID.append(MID); l_MID.append(MID);
                l_Sex.append(Sex); l_Sex.append(Sex);
                
            f_out.write(' '.join(l_FID)+'\n')
            f_out.write(' '.join(l_IID)+'\n')
            f_out.write(' '.join(l_PID)+'\n')
            f_out.write(' '.join(l_MID)+'\n')
            f_out.write(' '.join(l_Sex)+'\n')
            
            
    ### Main iteration
    count = 0

    while 1:
        
#         print(count)
        
        line_bgl = f_bgl.readline()
        
        if not bool(line_bgl):
            # EOF
            break
            
        else:
#             print(line_bgl)
#             print(line_bgl.startswith('M'))
            
            if not line_bgl.startswith('M'):
                # header line
                if not flag_useFAM:
                    f_out.write(line_bgl)
                    
            else:
                # marker line
                line_bim  = f_bim.readline()
#                 print(line_bim)
                [_chr, _label, _GD, _BP, _A1, _A2] = re.split(r'\s+', line_bim.rstrip('\n'))
                
                if not (_A1 in normal_nts and _A2 in normal_nts):
                    
                    """
                    T : _A1 (ALT)
                    A : _A2 (REF)
                    """
                                        
                    f_out.write(' '.join([(_A1 if item == 'T' else _A2 if item == 'A' else item) for item in re.split(r'\s+', line_bgl.rstrip('\n'))]) + '\n')
                    
                else:
                    f_out.write(line_bgl)
                    
                    
            count += 1
#             if count > 20 : break;

                    
            
    f_bgl.close()
    f_bim.close()
    f_out.close()
    
    return _out



if __name__ == '__main__':
    
    [_bgl, _ref_bim, _out, _tar_fam] = sys.argv[1:]