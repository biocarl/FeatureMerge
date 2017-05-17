from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Data.CodonTable import TranslationError
import sys, getopt, re, os, shutil, inspect
from exceptions import IndexError


#print to stderr (for piping)
def eprint(string) :
	print >> sys.stderr, string

#main
def main(argv):

	#Parameters
	input_file = ''
	input_file_corr = ''
	output_file = ''
        MAX_PAIR = 2 #Don't change!

	try:
		opts, args = getopt.getopt(argv,":ih",["input="])

		for opt, arg in opts:
			if opt == '-h':
				eprint ("See README")
				sys.exit()
	
			elif opt in ("-i","--input"):
				input_file = arg


		if (input_file == '') : 
			raise getopt.GetoptError ("Provide Files")

	except ( getopt.GetoptError, IndexError ) :
		eprint ('Wrong CLI-Arguments')
		sys.exit(2)
	

	#CODE
	###################################

        eprint(input_file)

        input_file_corr = str(os.path.splitext(input_file)[0])+'_corr.embl'
        output_file = str(os.path.splitext(input_file)[0])+'_merged.embl'


        #File-Handle
        input_file_corr_FH = open(input_file_corr , 'w')
        output_file_FH = open(output_file , 'w')
        input_file_FH = open(input_file , 'r')

        header = ["ID   DummyHeader; SV 1; linear; unassigned DNA; STD; UNC; 2572312 BP.","XX","AC   ---;","XX","DE   ---","XX","XX","OS   ---","OC   ---; --- ---; ----; ---;","OC   ---; ---; ---;","OC   ---","XX","FH   Key             Location/Qualifiers","FH"]

        #Transfer Dummy-Header
        input_file_corr_FH.write('\n'.join(header) + '\n')
        shutil.copyfileobj(input_file_FH, input_file_corr_FH)
        
        #Close FH
        input_file_corr_FH.close()
        input_file_FH.close()
        #output_file_FH.close()

        #Merge Features
        input_record = SeqIO.read(open(input_file_corr,"r"), "embl")
        features = input_record.features
        size = len(features)
        count_pair = 0 #<=MAX_PAIR
        is_uniq = True
        features_new = list()

        #LOG
        count_double = 0
        count_uniq = 0
        #/LOG

        for i in range(size):
            is_uniq=True
            count_pair = 1
            for j in range(i+1,size):
                if ((features[i].location.end == features[j].location.end)):
                    #Here is the pair j,i which is the same CDS : Create a new feature and merge the attributes
                    features[i].qualifiers.update(features[j].qualifiers)
                    features_new.append(features[i])
                    features[j].qualifiers.update({'duplicate':True})

                    #flags/counts
                    is_uniq = False
                    count_pair += 1 #count how many duplicates are allowed
                    count_double += 1
            if(is_uniq and not features[i].qualifiers.get('duplicate', False)):
                #Here is the unique CDS that will be transferred normally
                features_new.append(features[i])
                count_uniq += 1

            if(count_pair > MAX_PAIR):
                raise Exception("The file contains more than "+str(MAX_PAIR-1)+" duplicates of one feature")

        eprint(count_double)
        eprint(count_uniq)

        #Save merged features to the output file
        input_record.features = []

        input_record.features.extend(features_new)

        output_handle = open(output_file , 'w')
        SeqIO.write(input_record, output_handle, "embl")
        output_handle.close()

if __name__ == "__main__":

	main(sys.argv[1:])
