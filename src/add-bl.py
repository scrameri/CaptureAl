#!/usr/bin/env python
'''
Created on Jun 3, 2011
@author: smirarab
'''
import dendropy # needs dendropy version 3, version 4 will throw an error
import sys
import os
import copy
import os.path

if __name__ == '__main__':

    if len(sys.argv) < 3: 
        print "USAGE: [postfix|-|--] treefile"
        sys.exit(1)
    stdout = False
    if sys.argv[1] == "-":
        resultsFile = sys.stdout
        stdout = True
    elif sys.argv[1] == "--":
        postfix = "blen"
    else:
        postfix = sys.argv[1]
    
    c={}
    for treeName in sys.argv[2:]:
        if not stdout:
            resultsFile=open("%s.%s" % (treeName, postfix),'w')
        trees = dendropy.TreeList.get_from_path(treeName, 'newick', rooted=True) # S. Crameri, March 2019: removed rooted=True argument due to error message
        for tree in trees:
            for e in tree.postorder_edge_iter():
                if not e.length:
                    e.length = 1
        print >>sys.stderr, "writing results to " + resultsFile.name        
trees.write(resultsFile, 'newick', write_rooting=False) # S. Crameri, March 2019: removed 'newick' and write_rooting=False argument due to error message
