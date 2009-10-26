#!/u/tracyt/env/bin/python
import sys
# test 
import shutil
import cPickle
import traceback
import os, os.path
import time

import _mypath
import blastkit
import blastparser
from pygr import seqdb

import cgi

from bsddb import btopen
from shelve import BsdDbShelf

from operator import itemgetter

import fasta
import re                                 
import jinja2
###

# test modify~

PLACEHOLDER_MESSAGE = '''
<head><META http-equiv="refresh" content="5;"></head>
<body>
This page will automatically reload, and the results will be available here.
<p>
If you get a blank page, hit Reload on your browser.
Results for large datasets can take up to 20 minutes to be returned ...
</body>
'''

###

def do_cgi():
    """
    Main CGI function.  Retrieve form information, set up task, return
    placeholder, and spawn worker process.
    """
    
    # retrieve sequence from submitted form info
    form = cgi.FieldStorage()
   # name = form['name'].value
    sequence = form['sequence'].value
    
    input_filename = None
    if 'filename' in form:
        item = form['filename']
        input_filename = item.filename


    blastdb = form['blastdb'].value
    blast_program = form['blastprogram'].value
    e_value =float( form['e_value'].value)

    if "/" in blastdb:
        print "Content-type: text/html\n"
        print 'wrong argument!'
        sys.exit(2)




# check blastdb is a proper argument  --Q.Z.
 #   blastdb_outfile = blastdb + '.db'
    dbfile = '/u/qingpeng/blastkit/db/' + blastdb + '.fa'

    if not os.path.exists(dbfile):
        print "Content-type: text/html\n"
        print 'Database file does not exist. Please have a check!'
        sys.exit(2)


    # If neither a sequence or a file is submitted, request them

    if (not input_filename and not sequence):
        print "Content-type: text/html\n"
        print 'Please upload a file or submit a sequence'
        sys.exit(2)

    # If both are submitted, ask that they submit only one

    if (input_filename and sequence):
        print "Content-type: text/html\n"
        print 'Please only submit a file or a sequence'
        sys.exit(2)

    # make a working directory to save stuff in
    tempdir, dirstub = blastkit.make_dir()

    # If a single sequence is submitted / or multiple fasta sequences are submitted. -- Q.Z.
    
    if form.has_key("sequence"):
            #write out the query sequences with ">"
            fp = open('%s/query.fa' % (tempdir,), 'w')
            fp.write( '%s' %(sequence,))
            fp.close()
            

    # if it's a file that's uploaded instead
    if input_filename:
        data = item.file

        try:
            blast_input = open("%s/query.fa" % (tempdir), 'w')
        except:
            print "Content-type: text/html\n"
            print 'Cannot open', blast_input, 'for writing'
            sys.exit(2)


        try:
            fasta_data = fasta.load(data)
        except:
            print "Content-type: text/html\n"
            print 'This file does not seem to be a fasta file.  Please try again with a fasta file'
            sys.exit(0)

        for key in fasta_data:
            blast_input.write('>%s\n%s\n' % (key, fasta_data[key]))
        blast_input.close()




    # write out the placeholder message
    fp = open('%s/index.html' % (tempdir,), 'w')
    fp.write(PLACEHOLDER_MESSAGE)
    fp.close()

    # fork response function / worker function
    blastkit.split_execution(response_fn, (dirstub,), worker_fn, (tempdir,blastdb,blast_program,e_value))

def response_fn(dirstub):
    """
    Construct a new location & go there.
    """
    name = os.environ['SERVER_NAME']
    port = os.environ['SERVER_PORT']
    path = os.path.dirname(os.environ['SCRIPT_NAME'])
    
    url = 'http://%s:%s%s/files/%s/' % (name, port, path, dirstub)


    time.sleep(1)

    print 'Location: %s\n' % (url,)

            

@blastkit.write_tracebacks_to_file
def worker_fn(tempdir, blastdb, blast_program, e_value):
    """
    Run the BLAST and display the results.
    """

#    dbfile = '/u/tracyt/blastkit/db_raw/MA1W2.fa'
    dbfile = '/u/qingpeng/blastkit/db/' + blastdb + '.fa'
    newfile = tempdir + '/' + 'query.fa'
    str_evalue=str(e_value)
#  modify -b --Q.Z
    argu=["-e",str_evalue,"-b","500"]
    out, err = blastkit.run_blast(blast_program, newfile, dbfile,argu)

    fp = open(tempdir + '/blast-out.txt', 'w')
    fp.write(out)
    fp.close()

    if err and err.strip():
        fp = open(tempdir + '/blast-err.txt', 'w')
        fp.write(err)
        fp.close()

    results = list(blastparser.parse_string(out))

    fp = open(tempdir + '/blast-results.pickle', 'w')
    cPickle.dump(results, fp)
    fp.close()



    ###

    db = seqdb.SequenceFileDB(dbfile)
###########################

    display(tempdir,results,db,e_value,blastdb)

def display(tempdir,results,db,e_value,blastdb):
# a sub function to deal with the blast result and display



    hit_dict = {}
    hit_dict_query = {}
    hit_dict_blast = {}
    hit_dict_seq = {}

    
    fp = open(tempdir + '/index.html', 'w')
    for query in results:
        # For each match in the query
# if a subject matches to many queries, pick the match with lowest e-value

        for subject in query:
            flag = 0
            for hit in subject:
#                if (hit.expect < e_value):         
                    for keys in hit_dict:  # foreach subject that has been in the list, if has lower e-value, then don't add the newer record.                  
                        if subject.subject_name == keys and hit_dict[keys] < hit.expect:
                            flag = 1
                    if flag == 0:
                        hit_dict[subject.subject_name] = hit.expect
                        hit_dict_query[subject.subject_name] = query.query_name


        sorted_hits = sorted(hit_dict.iteritems(), key=itemgetter(1))

    
    of = open(tempdir + '/sequences.fa', 'w')
    xl = open(tempdir + '/xls.txt', 'w')
    xl.write('%s\t%s\t%s\t%s\n' %('Query_Seq','Target_Seq','e_value','Annotation_from_NR')) 
    blastdb_outfile = './db/'+blastdb + '.db'
    if os.path.exists(blastdb_outfile):
        file_test=True
        _db = btopen(blastdb_outfile)
        db2=BsdDbShelf(_db)
    else:
        file_test=False
        fp.write('Note: No annotation information available temporarily.<p>')    
#initiat 
    
    for unique_hits in sorted_hits:
            print_hit = unique_hits[0]
            seq = db[unique_hits[0]]
            seq = str(seq)

            of.write('>%s\n%s\n' % (print_hit, seq,))
            if file_test :
                    record = db2[print_hit]
                    if len(record) >0: # maybe 
                        for name,match_score in record:
                            xl.write('%s\t%s\t%s\t%s\n' %(hit_dict_query[print_hit], print_hit, unique_hits[1],name))
                            #splitPattern = re.compile('(.+?)\|(.+?)\|(.*?)(\[.*)')
                            #tmp_splitPattern  = re.compile("\[")
                            #if splitPattern.search(name):# if there is "[" in name 
                      
                             #   groups = splitPattern.search(name).groups()
    #                       the annotation may start with 'ref' or without 'ref'  
                              #  if groups[0] =='ref':
                               #     to_write = "ref|<a href=http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=retrieve&db=protein&list_uids="+groups[0]+"&dopt=full_report>"+groups[1]+"</a>      "+"<b>"+groups[2]+"</b>"+groups[3]
                                #else:
                         #           to_write=groups[0]+"|"+groups[1]+"<b>"+groups[2]+"</b>"+groups[3]  
                            #else: #if there is not any "[" in name
                          #      to_write = name 
                            
                            #fp.write('%s<br>' % (to_write,))
                    else:
                   # write to 'tab' delimiter file 
                        xl.write('%s\t%s\t%s\tN/A\n' %(hit_dict_query[print_hit], print_hit, unique_hits[1])) 
            else:
                xl.write('%s\t%s\t%s\tN/A\n' %(hit_dict_query[print_hit], print_hit, unique_hits[1]))



# template#################
# load template files from ./templates/
    thisdir = os.path.dirname(__file__)
    templatesdir = os.path.join(thisdir, 'templates')
    templatesdir = os.path.abspath(templatesdir)

    loader = jinja2.FileSystemLoader(templatesdir)
    env = jinja2.Environment(loader=loader)

                    # ok, now load the example template.
    template = env.get_template('result.html')


    html = template.render(locals())
    fp.write(html)
                                



    fp.close()
    of.close()
    xl.close()
###

try:
    do_cgi()
except SystemExit:
    pass
except:                                 # catch errors & write them out
    print 'Content-type: text/html\n'
    print '<pre>'
    print traceback.format_exc()
    print '</pre>'

