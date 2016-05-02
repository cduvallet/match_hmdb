import re

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='input "xml file", concatenated from many individual xml files', default='', required=True)
    parser.add_argument('-o', help='output xml file, with excess declarations removed', default='', required=True)
    args = parser.parse_args()
    inputxml = args.i
    outputxml = args.o
else:
    #Path to hmdb xml file
    inputxml = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/hmdb_metabolites/hmdb_metabolites.xml'
    #Path to your output file
    outputxml = '/home/irockafe/Documents/MIT/Fall_2015/Alm_Lab/METLIN_parse_hmdb_search/hmdb_metabolites/removed_xml_hmdb_metabolites.xml'

xml_top_tag = 'database'

with open(inputxml, 'r+') as f, open(outputxml, 'w') as output:
	#Add a new root to the xml - you can only have one root
	for line_num, line in enumerate(f):
		#If the line contains an xml header, don't write it to
		#the new file, except the first tag
		if re.search('<\?xml version', line) and line_num == 0:
			#write out the xml declaration and the first tag
			output.write(line)
			output.write('<'+xml_top_tag+'>\n')
		elif re.search('<\?xml version', line):
			pass
		else:
			#keep the indentation the same as the original file
			output.write('  '+line)
			#write out something for user to see it's working.
			if line_num % 1000000 == 0:
				print "Writing line %s" % line_num
	#Close the root
	output.write('</'+xml_top_tag+'>')
