import urllib2
import time

# do not use a file for running the script there was some error ! so use a list instead. with all ids in a single list
# you need to run the script with a variable Genelist 

Genelist=[#a list containing all the cpath id's#] 
 
for i in GeneList:
    result=open('%s.xml' %i,'a')
    response = urllib2.urlopen('http://www.pathwaycommons.org/pc/webservice.do?cmd=get_record_by_cpath_id&version=2.0&q=%s&output=biopax'%i)
    html= response.read()
    #print html
    result.write(html)
    result.close()
    time.sleep(1)




os.chdir('Network_analysis/Nicole_mTERF_und_Co/python_Id_conversion_out')