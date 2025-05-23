#!/usr/bin/env python
######################################################################
#make_kreport.py takes in the kraken output file and the make_ktaxonomy.py 
#output file to generate a kraken report file 
#Copyright (C) 2020 Jennifer Lu, jennifer.lu717@gmail.com
#
#This file is part of KrakenTools
#KrakenTools is free software; oyu can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.
#
######################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Updated: 04/15/2020
#
#This program creates the kraken report file from
#the make_ktaxonomy.py output and the kraken output file
#
#Required Parameters:
#   -i,-k,--kraken X....................kraken output file
#   -t,--taxonomy X.....................taxonomy file 
#   -o, --output X......................output kraken report file 
#Optional Parameters:
#   -h, --help..........................show help message.
#################################################################################
import os, sys, argparse
import operator
from time import gmtime
from time import strftime 
#################################################################################
#Tree Class
#usage: tree node used in constructing taxonomy tree
class Tree(object):
    'Tree node.'
    def __init__(self,  taxid, name, level_rank, level_num, p_taxid, parent=None,children=None):
        self.taxid = taxid
        self.name = name
        self.level_rank= level_rank
        self.level_num = int(level_num)
        self.p_taxid = p_taxid
        self.all_reads = 0
        self.lvl_reads = 0
        #Parent/children attributes
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)
#################################################################################
#Main method
def getArgs():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', '-k','--kraken', dest='kraken_file', required=True,
        help='Kraken output file (5 tab-delimited columns, taxid in 3rd column)')
    parser.add_argument('-t','--taxonomy', dest='tax_file', required=True,
        help='Output taxonomy file from make_ktaxonomy.py')
    parser.add_argument('-o','--output',dest='out_file', required=True,
        help='Output kraken report file')
    parser.add_argument('-c','--classifier', required=False, default='centrifuge',
        help='Classifier used, centrifuge, kraken2, or kraken2server [default: centrifuge]')
    parser.add_argument('--use-read-len',dest='use_read_len',
        action='store_true',default=False, required=False,
        help='Make report file using sum of read lengths [default: read counts]')
    args = parser.parse_args()
    return args

def readTaxFile(tax_file):
    #STEP 1/4: READ TAXONOMY FILE  
    count_nodes = 0
    #sys.stdout.write(">> STEP 1/4: Reading taxonomy %s...\n" % tax_file)
    #sys.stdout.write("\t%i nodes saved" % (count_nodes))
    #Parse taxonomy file 
    root_node = -1
    taxid2node = {}
    t_file = open(tax_file,'r')
    for line in t_file:
        count_nodes += 1
        #sys.stdout.write("\r\t%i nodes saved" % (count_nodes))
        #sys.stdout.flush()
        [taxid, p_tid, rank, lvl_num, name] = line.strip().split('\t|\t')
        curr_node = Tree(taxid, name, rank, lvl_num, p_tid)
        taxid2node[taxid] = curr_node
        #set parent/kids
        if taxid == "1":
            root_node = curr_node
        else:
            curr_node.parent = taxid2node[p_tid]
            taxid2node[p_tid].add_child(curr_node)
    t_file.close()
    sys.stdout.write("\r\t%i nodes saved\n" % (count_nodes))
    sys.stdout.flush()

    return taxid2node, root_node

def readKrakenFile(kraken_file, classifier='centrifuge',use_read_len=True):
    read_count = 0
    base_count = 0
    #sys.stdout.write(">> STEP 2/4: Reading kraken file %s...\n" % kraken_file)
    #sys.stdout.write("\t%i million reads processed" % read_count)
    #sys.stdout.flush()
    #Save counts per taxid
    taxid2counts = {}
    taxid2allcounts = {} 
    k_file = open(kraken_file,'r')
    for line in k_file:
        read_count += 1
        #if read_count % 1000 == 0:
            #sys.stdout.write('\r\t%0.3f million reads processed' % float(read_count/1000000.))
            #sys.stdout.flush()
        l_vals = line.strip().split('\t')
        if l_vals[0] == 'readID': continue
        if classifier=='centrifuge':
            base_count += int(l_vals[6])
        elif classifier in ['kraken2', 'kraken2server']:
            base_count = int(l_vals[3])
        
        taxid = l_vals[2]
        count = 1
        #If using read length instead of read counts
        if use_read_len:
            if '|' in l_vals[3]:
                [len1,len2] = l_vals[3].split('|')
                count = int(len1)+int(len2)
            else:
                if classifier=='centrifuge':
                    count = int(l_vals[6])
                elif classifier in ['kraken2', 'kraken2server']:
                    count = int(l_vals[3])
        #add to dictionaries 
        if taxid not in taxid2counts:
            taxid2counts[taxid] = count
            taxid2allcounts[taxid] = count
        else:
            taxid2counts[taxid] += count
            taxid2allcounts[taxid] += count
    k_file.close()
    #.stdout.write('\r\t%0.3f million reads processed\n' % float(read_count/1000000.))
    #sys.stdout.flush()
    return taxid2counts, taxid2allcounts, read_count, base_count

def calculateTaxIDreads(taxid2counts, taxid2node, taxid2allcounts):
    #sys.stdout.write(">> STEP 3/4: Creating final tree...\n")
    for curr_tid in taxid2counts:
        #Skip unclassified
        if curr_tid == '0':
            continue 
        #print(taxid2node)
        if curr_tid not in taxid2node: 
            print(f'{curr_tid} not in taxid2nodes')
            continue
        p_node = taxid2node[curr_tid].parent 
        add_counts = taxid2counts[curr_tid] 
        #Assign reads for node
        taxid2node[curr_tid].lvl_reads += add_counts 
        taxid2node[curr_tid].all_reads += add_counts 
        while (p_node != None):
            #Add child reads to parent node 
            p_taxid = p_node.taxid
            if p_taxid not in taxid2allcounts:
                taxid2allcounts[p_taxid] = add_counts
            else:
                taxid2allcounts[p_taxid] += add_counts
            p_node.all_reads += add_counts
            #Get next parent node
            p_node = p_node.parent
    return taxid2allcounts

def makeReport(taxid2counts,taxid2allcounts,out_file, read_count, base_count, root_node):
    #sys.stdout.write(">> STEP 4/4: Printing report file to %s...\n" % out_file)
    o_file = open(out_file,'w')
    #Write line for unclassified reads:
    if '0'  in taxid2counts:
        o_file.write("%6.2f\t" % (float(taxid2counts['0'])/float(base_count)*100))
        o_file.write("%i\t%i\t" % (taxid2counts['0'],taxid2counts['0']))
        o_file.write('U\t0\tunclassified\n')
    #Get remaining lines 
    parse_nodes = [root_node]
    while len(parse_nodes) > 0:
        curr_node = parse_nodes.pop(0)
        curr_tid = curr_node.taxid 
        #Print information for this level
        if curr_tid in taxid2allcounts:
            o_file.write("%6.2f\t" % (float(taxid2allcounts[curr_tid])/float(base_count)*100))
            o_file.write("%i\t" % taxid2allcounts[curr_tid])
        else:continue
        if curr_tid not in taxid2counts:
            o_file.write("0\t")
        else:
            o_file.write("%i\t" % taxid2counts[curr_tid])
        o_file.write("%s\t" % curr_node.level_rank)
        o_file.write("%s\t" % curr_tid)
        o_file.write(" "*curr_node.level_num*2 + curr_node.name + "\n")
        #Add children to list
        for child in sorted(curr_node.children, key=operator.attrgetter('all_reads')):
            if child.taxid not in taxid2allcounts:
                continue
            if taxid2allcounts[child.taxid] == 0:
                continue
            #Add to list 
            parse_nodes.insert(0,child)    
    o_file.close() 
    #End of program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    #sys.stdout.write("PROGRAM END TIME: " + time + '\n')
    

def main():
    args=getArgs()
    
    #Start Program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')

    #STEP 1/4: READ TAXONOMY FILE  
    taxid2node, root_node = readTaxFile(args.tax_file)

    #STEP 2/4: READ KRAKEN FILE FOR COUNTS PER TAXID
    taxid2counts, taxid2allcounts, read_count, base_count = readKrakenFile(args.kraken_file, classifier=args.classifier)
    
    #STEP 3/4: FOR EVERY TAXID PARSED, ADD UP TOTAL READS
    taxid2allcounts = calculateTaxIDreads(taxid2counts, taxid2node, taxid2allcounts)

    #STEP 4/4: PRINT REPORT FILE 
    makeReport(taxid2counts, taxid2allcounts, args.out_file, read_count, base_count, root_node)
    sys.exit(0)

#################################################################################
if __name__ == "__main__":
    main()
#################################################################################
##################################END OF PROGRAM#################################
#################################################################################

