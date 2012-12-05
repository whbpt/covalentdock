#!/usr/bin/python
import sys
import string



def findAtom(filename):
    """read in the mol2 file for ligands, find possible covalent bond builders.
    Return the atom names of the possible covalent bond builder."""

    f=open(filename)
    
    atom=False
    bond=False
    atomList=[]
    bondList=[]
    buff=[]
    
    for line in f:
        if '@<TRIPOS>' in line:
            if 'ATOM' in line:
                atom=True
                bond=False
                buff+=[line]
                continue
            elif 'BOND' in line:
                bond=True
                atom=False
                buff+=[line]
                continue
            else:
                bond=False
                atom=False
        if atom:
            atomList+=[line.split()]
        elif bond:
            bondList+=[line.split()]
        else:
            buff+=[line]
    f.close()

    connect=[[] for x in range(len(atomList)+1)]

    for line in bondList:
        ind=string.atoi(line[1])
        des=string.atoi(line[2])
        typ=line[3]
        connect[ind]+=[(des, typ)]
        connect[des]+=[(ind, typ)]

    atomType=[x[5].split('.')[0] for x in atomList]

    double=filter(lambda x: x[3]=='2', bondList)
    pair=[[string.atoi(x[1]),string.atoi(x[2])] for x in double]
    ctoc=filter(lambda x: ((atomType[x[0]-1]=='C') and (atomType[x[1]-1]=='C')) and (atomList[x[0]-1][5]!='C.ar') and (atomList[x[1]-1][5]!='C.ar'), pair)
    otoc=filter(lambda x: ((atomType[x[0]-1]=='C') and (atomType[x[1]-1]=='O')) or ((atomType[x[0]-1]=='O') and (atomType[x[1]-1]=='C')), pair)
#    print ctoc
#    print otoc
#initialize candidate set
    candidate=[]

#type 1
    for line in otoc:
        if atomType[line[0]-1]=='C':
            carb=line[0]
        else:
            carb=line[1]
        temp=filter(lambda x: (string.atoi(x[1])==carb) or (string.atoi(x[2])==carb), bondList)
        contact=[]
        for t in temp:
            if (string.atoi(t[1])==carb):
                contact+=[string.atoi(t[2])]
            else:
                contact+=[string.atoi(t[1])]
        for t in ctoc:
            if ((t[0] in contact) and (t[1] not in contact)):
                candidate+=[(1,t[1])]
            if ((t[0] not in contact) and (t[1] in contact)):
                candidate+=[(1,t[0])]

#type 2
    for line in otoc:
        if atomType[line[0]-1]=='C':
            carO=line[0]
        else:
            carO=line[1]
        ni=filter(lambda x: (x[1]=='1') and (atomType[x[0]-1]=='N'), connect[carO])
        ca=filter(lambda x: (x[1]=='1') and (atomType[x[0]-1]=='C'), connect[carO])
        if (len(ni)*len(ca)!=1):
            continue
        nitro=ni[0][0]
        carb2=ca[0][0]
        con_ni=[x[0] for x in connect[nitro]]
        con_ca2=[x[0] for x in connect[carb2]]
        carb3=[x for x in con_ni if (x in con_ca2) and (atomType[x-1]=='C') and (x!=carO)]

        if (len(carb3)==1):
            candidate+=[(2,carO)]


    #print "Candidate(s):"
    #print candidate

    name=filename.split('/')[-1].split('.')[0]

    outAtom=[]
    for i in range(1,len(atomList)+1):
        outAtom+=["%d %s %s %s %s %s\n" %(i,"%s%d" %(atomType[i-1],i),atomList[i-1][2],atomList[i-1][3],atomList[i-1][4],atomList[i-1][5])]
    oriBond=[' '.join(x)+"\n" for x in bondList]
    f_ori=file(name+'_unbond.mol2','w')
    for line in buff:
        f_ori.write(line)
        if '@<TRIPOS>ATOM' in line:
            for l in outAtom:
                f_ori.write(l)
        elif '@<TRIPOS>BOND' in line:
            for l in oriBond:
                f_ori.write(l)
            break
    f_ori.close()
    
    total=0
    fileList=[]
    
    for record in candidate:
        type=record[0]

        if type==1:
            index=record[1]
            mid=filter(lambda x: x[1]=='2', connect[index])[0][0]
            carb=filter(lambda y: atomType[y[0]-1]=='C',filter(lambda x: x[1]=='1', connect[mid]))[0][0]
            out=outAtom[index-1]
            #print 'ATOM '+out.split(' ')[1]
    

            carb_cor=[string.atof(x) for x in atomList[carb-1][2:5]]
            mid_cor=[string.atof(x) for x in atomList[mid-1][2:5]]
            index_cor=[string.atof(x) for x in atomList[index-1][2:5]]
            
            a=[]
            b=[]
            for i in range(0,3):
                a+=[index_cor[i]-mid_cor[i]]
                b+=[carb_cor[i]-mid_cor[i]]
            v=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
            length=(v[0]**2+v[1]**2+v[2]**2)**0.5
            v_pos=[x/length for x in v]

            bondBond=[]
            for line in bondList:
                if (string.atoi(line[1])==index and (string.atoi(line[2])==mid)) or (string.atoi(line[2])==index and (string.atoi(line[1])==mid)):
                    bondBond+=[[line[0],line[1],line[2],'1']]
                else:
                    bondBond+=[line]

            bondAtom=[]
            for i in range(1,len(atomList)+1):
                if (i==mid) or (i==index):
                    bondAtom+=["%d %s %s %s %s %s\n" %(i,"%s%d" %(atomType[i-1],i),atomList[i-1][2],atomList[i-1][3],atomList[i-1][4],'C.3')]
                else:
                    if (atomList[i-1][1]=='DEL') or (atomList[i-1][1]=='Had') or (atomList[i-1][1]=='Dum'):
                        bondAtom+=["%d %s %s %s %s %s\n" %(i,atomList[i-1][1],atomList[i-1][2],atomList[i-1][3],atomList[i-1][4],atomList[i-1][5])]
                    else:
                        bondAtom+=["%d %s %s %s %s %s\n" %(i,"%s%d" %(atomType[i-1],i),atomList[i-1][2],atomList[i-1][3],atomList[i-1][4],atomList[i-1][5])]

            add_mid=[]
            pos=[]
            neg=[]
            for i in range(0,3):
                pos+=[v_pos[i]*1.09+mid_cor[i]]
                neg+=[v_pos[i]*-1.09+mid_cor[i]]
            add_mid+=['%.4f %.4f %.4f H' %(pos[0],pos[1],pos[2])]
            R1=filter(lambda x: (x[0]!=index) and (x[0]!=carb), connect[mid])
            if (len(R1)==0 ) or (atomType[R1[0][0]-1]!='H'):
                add_mid+=['%.4f %.4f %.4f H' %(neg[0],neg[1],neg[2])]

#    print add_mid

            add_index=[]
            pos=[]
            neg=[]
            for i in range(0,3):
                pos+=[v_pos[i]*1.83+index_cor[i]]
                neg+=[v_pos[i]*-1.83+index_cor[i]]
            add_index+=['%.4f %.4f %.4f S:%.4f %.4f %.4f H' %(pos[0],pos[1],pos[2],pos[0]+v_pos[0]*1.31, pos[1]+v_pos[1]*1.31, pos[2]+v_pos[2]*1.31)]
            R23=filter(lambda x: x[0]!=mid, connect[index])
            if (atomType[R23[0][0]-1]!=atomType[R23[1][0]-1]):
                add_index+=['%.4f %.4f %.4f S:%.4f %.4f %.4f H' %(neg[0],neg[1],neg[2],neg[0]-v_pos[0]*1.31, neg[1]-v_pos[1]*1.31, neg[2]-v_pos[2]*1.31)]
#    print add_index
            total=0
            atom_count=len(bondAtom)+3
            bond_count=len(bondBond)+3
            for aH in add_mid:
                for aS in add_index:
                    total+=1
                    f_now=file(name+"_bond_1_%s_%d.mol2" %(out.split(' ')[1],total),"w")
                    fileList+=[name+"_bond_1_%s_%d.mol2" %(out.split(' ')[1],total)]
        
                    f_now.write('@<TRIPOS>MOLECULE\n')
                    f_now.write(name+"_bond_1_%s_%d\n" %(out.split(' ')[1],total))
                    f_now.write("%d %d 1 0 0\n" %(atom_count,bond_count))
                    f_now.write("SMALL\n")
                    f_now.write("USER_CHARGES\n\n")
                    f_now.write("@<TRIPOS>ATOM\n")
                    for l in bondAtom:
                        f_now.write(l)
                    f_now.write("%d %s %s\n" %(len(bondAtom)+1, "Had" ,aH))
                    f_now.write("%d %s %s\n" %(len(bondAtom)+2, "Dum" ,aS.split(':')[0]))
                    f_now.write("%d %s %s\n" %(len(bondAtom)+3, "DEL" ,aS.split(':')[1]))
                    f_now.write("@<TRIPOS>BOND\n")
                    for l in bondBond:
                        f_now.write(' '.join(l)+"\n")
                    f_now.write("%d %d %d 1\n" %(len(bondBond)+1, mid, len(bondAtom)+1))
                    f_now.write("%d %d %d 1\n" %(len(bondBond)+2, index, len(bondAtom)+2))
                    f_now.write("%d %d %d 1\n" %(len(bondAtom)+3,len(bondAtom)+2,len(bondAtom)+3))
                    f_now.close()
        elif type==2:
            carb1=record[1]
            out=outAtom[carb1-1]
            nitro=filter(lambda x: (x[1]=='1') and (atomType[x[0]-1]=='N'), connect[carb1])[0][0]
            carb2=filter(lambda x: (x[1]=='1') and (atomType[x[0]-1]=='C'), connect[carb1])[0][0]
            con_ni=[x[0] for x in connect[nitro]]
            con_ca2=[x[0] for x in connect[carb2]]
            carb3=[x for x in con_ni if (x in con_ca2) and (atomType[x-1]=='C') and (x!=carb1)][0]

            nitro_cor=[string.atof(x) for x in atomList[nitro-1][2:5]]
            carb1_cor=[string.atof(x) for x in atomList[carb1-1][2:5]]
            carb2_cor=[string.atof(x) for x in atomList[carb2-1][2:5]]
            carb3_cor=[string.atof(x) for x in atomList[carb3-1][2:5]]

#add H-N
            v1=[]
            for i in range(0,3):
                v1+=[nitro_cor[i]-carb3_cor[i]]
            v1_mod=(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2])**0.5
            v1_std=[x/v1_mod for x in v1]
            
            add_H=[]
            for i in range(0,3):
                add_H+=[nitro_cor[i]+v1_std[i]*1.0]

#add O=C-OH
            v2=[]
            for i in range(0,3):
                v2+=[carb1_cor[i]-carb2_cor[i]]
            v2_mod=(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2])**0.5
            v2_std=[x/v2_mod for x in v2]
            
            add_O=[]
            add_OH=[]
            for i in range(0,3):
                add_O+=[carb1_cor[i]+v2_std[i]*1.4]
                add_OH+=[carb1_cor[i]+v2_std[i]*(1.4+1.0)]

            bondAtom=[]
            for i in range(1,len(atomList)+1):
                if (atomList[i-1][1]=='DEL') or (atomList[i-1][1]=='Had') or (atomList[i-1][1]=='Dum'):
                    bondAtom+=["%d %s %s %s %s %s\n" %(i,atomList[i-1][1],atomList[i-1][2],atomList[i-1][3],atomList[i-1][4],atomList[i-1][5])]
                else:
                    bondAtom+=["%d %s %s %s %s %s\n" %(i,"%s%d" %(atomType[i-1],i),atomList[i-1][2],atomList[i-1][3],atomList[i-1][4],atomList[i-1][5])]

            bondBond=[]
            bond_count=0
            for line in bondList:
                if (string.atoi(line[1])==carb1 and (string.atoi(line[2])==nitro)) or (string.atoi(line[2])==carb1 and (string.atoi(line[1])==nitro)):
                    continue
                else:
                    bondBond+=["%d %s %s %s\n" %(bond_count+1,line[1],line[2],line[3])]
                    bond_count+=+1

            total=0
            atom_count=len(atomList)
            total+=1
            f_now=file(name+"_bond_2_%s_%d.mol2" %(out.split(' ')[1],total), "w")
            fileList+=[name+"_bond_2_%s_%d.mol2" %(out.split(' ')[1],total)]

            f_now.write('@<TRIPOS>MOLECULE\n')
            f_now.write(name+"_bond_2_%s_%d.mol2\n" %(out.split(' ')[1],total))
            f_now.write("%d %d 1 0 0\n" %(atom_count+3,bond_count+3))
            f_now.write("SMALL\n")
            f_now.write("USER_CHARGES\n\n")
            f_now.write("@<TRIPOS>ATOM\n")
            for l in bondAtom:
                f_now.write(l)
            f_now.write("%d %s %.4f %.4f %.4f H\n" %(atom_count+1, "Had" , add_H[0],add_H[1],add_H[2]))
            f_now.write("%d %s %.4f %.4f %.4f O\n" %(atom_count+2, "Dum" , add_O[0],add_O[1],add_O[2]))
            f_now.write("%d %s %.4f %.4f %.4f H\n" %(atom_count+3, "DEL" , add_OH[0],add_OH[1],add_OH[2]))

            f_now.write("@<TRIPOS>BOND\n")
            for l in bondBond:
                f_now.write(l)
            f_now.write("%d %d %d 1\n" %(bond_count+1, nitro, atom_count+1))
            f_now.write("%d %d %d 1\n" %(bond_count+2, carb1, atom_count+2))
            f_now.write("%d %d %d 1\n" %(bond_count+3, atom_count+2, atom_count+3))
            f_now.close
            
    for line in fileList:
                print line

    return fileList

if __name__=="__main__":
    filename=sys.argv[1]
#    print "working on "+filename+"."
#    print filename
    findAtom(filename)

