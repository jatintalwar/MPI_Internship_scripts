############# script to edit the output file from parsing all the xml files in a more readable format   ############

import os 
os.chdir('/home/talwar/Network_analysis/Nicole_mTERF_und_Co/xml_out_analysis/')


new=open('output_edited.txt','a')
F=open('output_new_format.txt','r')
for line in F:
    Id=['NA']*4                                      # initialise list for Id values for every line in the file 
    Db=['NA']*4                                      # initialise list for DataBase values for every line in the file
    leftover=line.split(' ')                         # first split for spaces between lines for each individual entry (one line in output file was a unique entry)
    for i in leftover:
        if "ENTREZ_GENE" in i:						 # now i search for the presence of these strings in the individual elements after initial split command
            Db[0],Id[0]=i.split(':')                
        elif "GENE_SYMBOL" in i:
            Db[1],Id[1]=i.split(':')
        elif "CPATH" in i:							 # and if present the string is further split and the ID values are stored at a particular index in the list 
            Db[2],Id[2]=i.split(':')
        elif "UNIPROT" in i:
            Db[3],Id[3]=i.split(':')
    new.write("%s\n" % Id) #new.write("\n")          # later i print all the values inside the list "ID" and as i printed them at a particular index i know which value 
  													 # corresponds to which ID because they are all printed in the same way every time.

new.close()

###### there were multiple 'UNIPROT' & "ENTREZ_GENE" in the list for a single CPATH ID and GENE_SYMBOL hence i took only one of those through this script. 
###### i am guessing that is the first one. 

###### script was checked on 8/10/14 and was working fine like this !
