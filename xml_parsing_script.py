##############   code for parsing an xml file only one file here .. look down for the code for multiple files:   #################


from xml.dom import minidom

xmldoc=minidom.parse('72860.xml')

RDF = xmldoc.getElementsByTagName('rdf:RDF')[0]     # [0] is somehow important to put : doesnt work otherwise 
protein= RDF.getElementsByTagName("bp:protein")[0]

########   we need all bp:XREF tags in the xml file: it contains info about uniprot,ENTREZ_ID,Gene_symbol etc ############

xrefs= protein.getElementsByTagName("bp:XREF")

result_file=open('output.txt','a')
for xref in xrefs:
    relationship= xref.getElementsByTagName("bp:relationshipXref")
    for r in relationship:
        db=r.getElementsByTagName("bp:DB")[0].firstChild.data
        id=r.getElementsByTagName("bp:ID")[0].firstChild.data
        print(db + ":" + id)
        result_file.write(db+";"+id)
        result_file.write('\n')

for xref in xrefs:
    unification= xref.getElementsByTagName("bp:unificationXref")
    for u in unification:
        db_uni=u.getElementsByTagName("bp:DB")[0].firstChild.data
        id_uni=u.getElementsByTagName("bp:ID")[0].firstChild.data
        result_file.write(db_uni + ":" + id_uni)
        result_file.write("\n")


result_file.write("#####-----------######")
result_file.close()



############## the script above works perfectly now for a particular file!! do not make any new changes! ####################

############## script for all files:   ##################################

for file in os.listdir("/home/talwar/Network_analysis/Nicole_mTERF_und_Co/xml_out_analysis"):
    xmldoc=minidom.parse(file)
    RDF = xmldoc.getElementsByTagName('rdf:RDF')[0]     # [0] is somehow important to put : doesnt work otherwise 
    protein= RDF.getElementsByTagName("bp:protein")[0]
    xrefs= protein.getElementsByTagName("bp:XREF")
    result_file=open('output.txt','a')
    
    for xref in xrefs:
        relationship= xref.getElementsByTagName("bp:relationshipXref")
        for r in relationship:
            db=r.getElementsByTagName("bp:DB")[0].firstChild.data
            id=r.getElementsByTagName("bp:ID")[0].firstChild.data
            #print(db + ":" + id)
            result_file.write(db+";"+id)
            result_file.write('\n')
    for xref in xrefs:
        unification= xref.getElementsByTagName("bp:unificationXref")
        for u in unification:
            db_uni=u.getElementsByTagName("bp:DB")[0].firstChild.data
            id_uni=u.getElementsByTagName("bp:ID")[0].firstChild.data
            #print(db_uni + ":" + id_uni)
            result_file.write(db_uni + ":" + id_uni)
            result_file.write("\n")
            
            
    result_file.write("#####-----------######")
    result_file.write('\n')



result_file.close()             ### very important to close the file !!! 
        
#################### the script works now for 191 files or so ! yahoooooo...............##########################
    
    
    