from cov_radii import *

def file_into_list(file):
    '''Takes a file and converts each line into an item of a list'''
    fatty_list=[]
    with open (file, 'r') as fatty_file:
        for line in fatty_file:
            fatty_list.append(line[:-1])

    for i in range(len(fatty_list),0,-1):
        if '.' not in fatty_list[i-1] or ('C' not in fatty_list[i-1] and 'H' not in fatty_list[i-1] and 'O' not in fatty_list[i-1]):
            fatty_list.pop(i-1)

    return(fatty_list)

def parse_line(number, line_of_file):
    '''Takes an integer and line of xyz file and returns that integer, chemical symbol, and list of atomic coordinates'''
    subitem=line_of_file.split()
    for i in range(1,len(subitem)):
        subitem[i]=float(subitem[i])
    return([number,subitem[0],(subitem[1:])])

def distance(list1, list2):
    '''takes two lists of coordinates of two atoms and returns the distance between the two'''
    difference=[list1[0]-list2[0],list1[1]-list2[1],list1[2]-list2[2]]
    distance=(difference[0]**2+difference[1]**2+difference[2]**2)**0.5
    return(distance)

def coordination(list1,list2):
    '''takes an atom's coordinates and returns a list of atoms that are bonded to the atom - atom is bonded if interatomic distance (calculated using function distance()) is less than sum of atomic radii + 0.4 A'''
    coordination_list=[]
    distance_list=[]
    for i in range(len(list2)):
        if list1!=list2[i]:
            interatom_dist=distance(list1[2],list2[i][2])
            interatom_dist_limit=cov_dictionary[list1[1]]+cov_dictionary[list2[i][1]]+0.4
            if interatom_dist_limit>=interatom_dist:
                coordination_list+=[list2[i]]
    return(coordination_list)

def is_acid(list1,list2):
    '''Takes a C atom and returns whether that C is part of a carboxyl functionality - checks whether the C atom is bonded to another C and two O'''
    criterion_c=False
    criterion_o=0
    acid=False
    if 'C' in list1[1] and len(list2)>=3:
            for item in list2:
                if 'C' in item[1]:
                    criterion_c=True
                elif 'O' in item[1]:
                    criterion_o+=1
            if criterion_c==True and criterion_o==2:
                acid=True
    else:
        return("Entered atom isn't a carbon atom")
        
    return(acid)

def is_terminal(list1,list2):
    '''Takes a C atom and returns whether that C is terminal in the aliphatic chain - checks whether the atom is bonded only to one C (no O)'''
    terminal=False
    c_count=0
    o_count=0
    if 'C' in list1[1]:
        for item in list2:
            if 'C' in item[1]:
                c_count+=1
            elif 'O' in item[1]:
                o_count+=1
        if c_count==1 and o_count==0:
            terminal=True
    else:
        return("Entered atom isn't a carbon atom")
        
    return(terminal)

def unsaturation(list1, list2):
    '''Takes a C atom and returns saturation level of atom - 0 if saturated, 1 if involved in one double bond, or 2 if involved in two double bonds or a triple bond'''
    if 'C' in list1[1]:
        no_bonds=len(list2)
        unsaturation_no=4-len(list2)
        return(unsaturation_no)
    else:
        return ("Entered atom isn't a carbon atom")


if __name__=="__main__":

    print("xyz to SMILES converter \n------------------------------------------")

    file_name=input("Please enter your file name with .xyz at the end. Make sure your file is in the same folder as this code or enter the path to the file. \n")
    fatty_file=file_into_list(file_name)
    
    parse_line_list=[]
    distance_list=[]
    
    for i in range(len(fatty_file)):
        parse_line_list+=[parse_line(i+1,fatty_file[i])]
    
    for i in range(len(parse_line_list)-1):
        distance_list+=[[parse_line_list[i][1],parse_line_list[i+1][1],distance(parse_line_list[i][2],parse_line_list[i+1][2])]]
    
    smiles_atoms=[]
    for i in parse_line_list:
        if 'H' not in i:
            smiles_atoms+=[i]
    
    #finds the C atom in carboxylic acid functionality
    for i in range(len(smiles_atoms)):
        check_acid=is_acid(smiles_atoms[i],coordination(smiles_atoms[i],parse_line_list))
        check_terminal=is_terminal(smiles_atoms[i],coordination(smiles_atoms[i],parse_line_list))
        #Once the acidic C is found, determines which O is the terminal and which is hydroxyl O
        if check_acid==True:
            acid_C=smiles_atoms[i]
            acid_bonds=coordination(smiles_atoms[i],parse_line_list)
            for j in acid_bonds:
                if 'O' in j:
                    maybe_OH=(coordination(j,parse_line_list))
                    if len(maybe_OH)==1:
                        CO_O=j
                    elif len(maybe_OH)==2:
                        OH_O=j
    
        elif check_terminal==True:
            terminal_C=smiles_atoms[i]
    
    smiles_list=[]
    
    #smiles_list populated with characteristic atoms that have already been identified (OH oxygen, C of C-OH, C of C=O, and terminal C) in correct order
    smiles_list+=[OH_O]+[acid_C]+[CO_O]+[terminal_C]  
    
    #smiles_atom is all non-H atoms in molecule
    
    item_number=len(smiles_list)-1
    
    for k in range(len(smiles_atoms)):
        added=False
        for i in range(item_number,-1,-1):
            if 'C' in smiles_list[i] and added==False:
                for j in coordination(smiles_list[i], smiles_atoms):
                    if 'C' in j and j not in smiles_list:
                        smiles_list[i:i]=[j]
                        added=True
                        item_number=i
    
    saturation_list=[]
    for i in smiles_list:
        if type(unsaturation(i,coordination(i,parse_line_list)))==int:
            saturation_list+=[unsaturation(i,coordination(i,parse_line_list))]
    
    smiles_bonds=[]
    
    for i in smiles_list:
        if 'O' not in i:
            smiles_bonds+=[i]
    
    counter=0
    for i in range(len(saturation_list)):
        
        if saturation_list[i]==1 and saturation_list[i-1]==1:
            smiles_bonds[i+counter:i+counter]=['=']
            counter+=1
        elif saturation_list[i]==2 and saturation_list[i-1]==1:
            smiles_bonds[i+counter:i+counter]=['=']
            counter+=1
        
        elif i<len(saturation_list)-1:
            if saturation_list[i]==2 and saturation_list[i-1]==2 and saturation_list[i+1]==0:
                smiles_bonds[i+counter:i+counter]=['#']
                counter+=1
    
            elif saturation_list[i]==2 and saturation_list[i-1]==2 and saturation_list[i+1]==1:
                smiles_bonds[i+counter:i+counter]=['=']
                counter+=1
        
        elif i==len(saturation_list)-1:
            if saturation_list[i]==2 and saturation_list[i-1]==2:
                smiles_bonds[i+counter:i+counter]=['#']
                counter+=1
    
    smiles_bonds[0:0]=[OH_O]
    smiles_bonds[2:2]=[CO_O]
    
    smiles=''
    for i in range(len(smiles_bonds)):
        if smiles_bonds[i]==CO_O:
            smiles+='(='+smiles_bonds[i][1]+')'
        elif len(smiles_bonds[i])==1:
            smiles+=smiles_bonds[i]
        else:
            smiles+=smiles_bonds[i][1]
    
    print("\n", smiles, "\n")
    
    input('Click enter to close.')
