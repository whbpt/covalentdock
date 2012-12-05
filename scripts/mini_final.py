#!/usr/bin/python
import sys
import string

def correctCord(basename):
    """read in the mol2 file for ligands, find possible covalent bond builders.
    Return the atom names of the possible covalent bond builder."""

#    print basename

    min_f=open(basename+'_amber.mol2')
    ori_f=open(basename+'.mol2')

    buff=[]

    ori_atom=[]
    min_atom=[]
    bondList=[]

    atom=False
    bond=False

    for line in min_f:
        if '@<TRIPOS>' in line:
            if 'ATOM' in line:
                atom=True
                bond=False
                continue
            elif 'BOND' in line:
                bond=True
                atom=False
                continue
            else:
                bond=False
                atom=False
        if atom:
            min_atom+=[line.split()]
        elif bond:
            bondList+=[line.split()]
        else:
            continue
    min_f.close()

    atom=False
    bond=False
    for line in ori_f:
        if '@<TRIPOS>' in line:
            if 'ATOM' in line:
                atom=True
                bond=False
                continue
            elif 'BOND' in line:
                bond=True
                atom=False
                continue
            else:
                bond=False
                atom=False
        if atom:
            ori_atom+=[line.split()]
        elif bond:
            continue
        else:
            continue
    ori_f.close()

    out_atom=[]
    out_atom_index=[]
    for i in range(0,len(ori_atom)):
        if (ori_atom[i][1]!='DEL'):
            out_atom+=[[ori_atom[i][0],ori_atom[i][1],min_atom[i][2],min_atom[i][3],min_atom[i][4],ori_atom[i][5]]]
            out_atom_index+=[ori_atom[i][0]]
        else:
            continue
    lut=[]
    for i in range(0,len(ori_atom)):
        if ori_atom[i][1]!='DEL':
            lut+=[out_atom_index.index(ori_atom[i][0])+1]
        else:
            lut+=[0]
#    for line in out_atom:
#        print line
    DEL=[x[0] for x in ori_atom if x[1]=='DEL']
    Had=filter(lambda x: x[1]=='Had', out_atom)
    Had_index=[x[0] for x in Had]
    #lut=[(x[0],"%d" %(out_atom.index(x)+1)) for x in out_atom]
    
    out_bond=[x for x in bondList if ((x[1] not in DEL) and (x[2] not in DEL))]
#    for line in out_bond:
#        print line

    f_out=open(basename+'_mini.mol2','w')

    f_out.write('@<TRIPOS>MOLECULE\n')
    f_out.write(basename+'\n')
    f_out.write('%d %d 1 0 0\n' %(len(out_atom), len(out_bond)))
    f_out.write('SMALL\n')
    f_out.write('gas\n\n')
    f_out.write('@<TRIPOS>ATOM\n')
    for i in range(0,len(out_atom)):
        f_out.write('%d %s\n' %(i+1,' '.join(out_atom[i][1:])))
    f_out.write('@<TRIPOS>BOND\n')
    for i in range(0,len(out_bond)):
        f_out.write("%d %d %d %s\n" %(i+1, lut[string.atoi(out_bond[i][1])-1],lut[string.atoi(out_bond[i][2])-1],out_bond[i][3]))
    f_out.close()

if __name__=="__main__":
    basename=sys.argv[1]
#    print "working on "+filename+"."
#    print filename
    correctCord(basename)

