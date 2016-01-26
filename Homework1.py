#"""This assignment is designed to introduce the genomic file formats and provide
#an introduction to the requirements for programming assignments.  The
#application will read a genome file in FASTA format, calculate the GC content
#across an entire chromosome, search the chromosome for exact matches to a given
#k-mer and report those matches using an annotation file format.  The k-mer
#match positions should be written to standard output in a BED or GFF format.
#
#The genome file and k-mer must be specified on the command line using the
#following format:
#    -f <filename> -c <chromo name> -k <k-mer>.
#
#    Requirements:
#        Process command line arguments to specify the genome file, the
#        chromosome to process, and the k-mer to match
#        Read genome file
#        Calculate %GC across specified chromosome
#        Search the genome for a specific k-mer (use sequence ATG as sample)
#        Report match locations using an annotation file format"""

import getopt, sys, os


def Usage():
    print "\nApplication: %s\n%s [options] -c <file> <sequences>      \n"    \
            "     -f <file1>    - File to be Read           \n"    \
            "     -c            - Chromosome to be read           \n"    \
            "     -k            - Kmer to search for in the Chromosme File   \n"    \
            "\n" %(sys.argv[0],sys.argv[0])

def Main():
        fflag = None
        cflag = None
        kflag = None
        try:    
            # getopt returns a list of options and a list of other args.
	    # Options listed in the 2nd parameter that require an argument
            #    are followed by ':'; the third parameter lists the multiple character
	    # 	 options such as "--help".
            options, other_args = getopt.getopt(sys.argv[1:], "f:c:k:", ["help"])
	
        except getopt.GetoptError:
            print "There were errors in parsing the command line!"
	    Usage()
            sys.exit(1)

	# process the list of options returned by getopt
	for option, value in options:
            if option == '-f':			
                fflag = value
	    elif option == '-c':
                cflag = value
	    elif option == '-k':
                kflag = value
	    elif option == '--help':
                Usage()
	    else: 
                print "Unhandled opt [%s][%s]"%(option)

	if (fflag == None):
            print "\n HEY! A file must be chosen. Parameter required."
	    Usage()
	    exit(-1)
		
	# place the processing here
        def GenomeReader(GenomeFile):
            """
            This functions purpose is to read in a header delineated sequence
            file. The output to this function is a dictionary that has the Key
            being the Header of the sequence and value as all line following
            the header until the next header symbol (>) us found
            """
            key = []
            GenomeScaffolds = {}
            with open(GenomeFile, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        NamedSeq = line.replace('>', '')
                        key.append(NamedSeq)
                        GenomeScaffolds[NamedSeq] = ""
                    else:
                        GenomeScaffolds[NamedSeq] += line
            return GenomeScaffolds
        

        def GCcounter(ChromosomeName, GenomicReadFile):
            """
            Given a specified chromsome from the genome file this file parses
            look in the dictiionary created above for the appropriate key and
            then analyzes the nucleortide strings thin the value of the
            dictionary. 
            """
            for key, value in GenomicReadFile.iteritems():
                if key == ChromosomeName:
                    TotalLen = len(value)
                    Counter = 0
                    for item in value:
                        if item == "G" or item == "C":
                            Counter += 1
                    PercentageGC = float(Counter)/float(TotalLen)
                    return PercentageGC * 100
                        
        
        def KmerSearcher(KMERToSearch, GenomicReadFile):
            """
            This function searches all the nucleotide sequence in a given
            dictionary for a KMER. It uses a sliding window approach to search
            every position of the genome for the correct kmer while the counter
            used to record location is less than the length of the scaffold. 
            If found the kmer location is returned as a list with a + or - to
            identify which strand the kmer lies on
            """
            KMERToSearch = str(KMERToSearch)
            SeqRevComp = RevComp(KMERToSearch)
            DNALOCATIONSTORAGE = []
            LengthOfKmer = len(KMERToSearch)
            for key, value in GenomicReadFile.iteritems():
                Counter2 = 0
                ScaffoldLen = int(len(value))
                while Counter2 < ScaffoldLen:
                    ChromosomeList = []
                    ChromosomeList.append(key)
                    for item in value:
                        if str(KMERToSearch) == str(value[Counter2:Counter2 + LengthOfKmer]):
                            ChromosomeList.append([Counter2, '+'])
                            Counter2 += 1
                        elif SeqRevComp == str(value[Counter2:Counter2 +
                            LengthOfKmer]):
                            ChromosomeList.append([Counter2, '-'])
                            Counter2 += 1 
                        else: 
                            Counter2 += 1
                    DNALOCATIONSTORAGE.append(ChromosomeList)
            return DNALOCATIONSTORAGE


        def RevComp(DNASTRING):
            """Reverse complements a DNA string in order to seach for both +
            and - direction of the given kmer. Will return a string with the
            reverse complement
            """
            NewRevComp = []
            DNASTRING = DNASTRING[::-1]
            for item in DNASTRING:
                if item == 'A':
                    NewRevComp.append('T')
                elif item == 'T':
                    NewRevComp.append('A')
                elif item == 'C':
                    NewRevComp.append('G')
                elif item == 'G':
                    NewRevComp.append('C')
            return ''.join(NewRevComp)

        def GBBWriter(LocationalityData, KMER):
            """ Writes a GFF file with the name OutputFile.gbb. Writes
            in standard format with all the above location data used.
            """
            KmerLengthToAdd = len(KMER)
            with open("OutputFile.gbb", 'a') as f:
                for item in LocationalityData:
                    ChromName = item[0]
                    ProgramName = "Pablo"
                    KmerName = str(KMER)
                    for thing in item:
                        if type(thing[0]) != int:
                            pass
                        else:
                            SeqStart = int(thing[0]) + 1
                            SeqEnd = int(thing[0]) + KmerLengthToAdd
                            Strand = thing[1]
                            RandomSeq = "ID=Name=%s;Pstart=%sPend=%s" %  (KmerName, SeqStart, SeqEnd)
                            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                                    ChromName, ProgramName, KmerName, SeqStart,
                                    SeqEnd, "100", Strand, '.', RandomSeq))

""" Below is execution of the code preceding"""

        THING = GenomeReader(fflag)
        PercentTEST = GCcounter(cflag, THING)
        print " \n The GC percentage over the Given chromosome %s was %s" % (cflag,
                PercentTEST)
        Locations = KmerSearcher(kflag, THING)
        GBBWriter(Locations, kflag)




	
	# return success or failure of the processing
	
###############################################################################
#  The following code will call the Main() routine if this code was called from 
#  the command line.  If this code was imported into another application, the
#  main function will not be called.
###############################################################################
if __name__ == '__main__':
	Main()
	
