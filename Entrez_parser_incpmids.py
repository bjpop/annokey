#!/bin/env python
import sys
#from lxml import etree
from Bio import Entrez
import re

def main():
    Hsxml = open('Homo_sapiens.211112xml')
    records = Entrez.parse(Hsxml)
    for record in records:
        status = record['Entrezgene_track-info']['Gene-track']['Gene-track_status']
        if status.attributes['value']=='discontinued':
            continue
        else:
            geneid = GeneID(record)
            genename = GeneName(record)
            synonyms = Syns(record)
            name = Name(record)
            altnames = AltNames(record)
            summary = Summary(record)
            generifs = GeneRIFs(record)
            function = Function(record)
            component = Component(record)
            process = Process(record)
            pathway = Pathway(record)
            interaction = Interaction(record)
            conserved = Conserved(record)
            pmids = PMID(record)
            outputlist = [geneid,genename,synonyms,name,altnames,summary,generifs,function,component,process,pathway,interaction,conserved,pmids]
            print "\t".join(str(thing) for thing in outputlist)

def GeneID(record):
    try:
        geneid = record['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
        return geneid
    except:
        return None
def GeneName(record):
    try:
        genename = record['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
        return genename
    except:
        return None

def Syns(record):
    try:
        synonyms = record['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
        return synonyms 
    except:
        return None

def Name(record):
    try:
        name = record['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
        return name
    except:
        return None

def AltNames(record):
    try:
        altnames = record['Entrezgene_prot']['Prot-ref']['Prot-ref_name']
        return altnames
    except:
        return None

def Summary(record):
    try:
        summary = record['Entrezgene_summary']
        return summary
    except:
        return None

def GeneRIFs(record):
    generiflist = []
    try:
        generifs = record['Entrezgene_comments']
        for item in generifs:
            for key, rif in item.iteritems():
                if key == 'Gene-commentary_text':
                    generiflist.append(rif)
        return generiflist
    except:
        return None

def Pathway(record):
    pathwaylist = []
    try:
        pathway = record['Entrezgene_comments']
        for item in pathway:
            for key,value in item.iteritems():
                if value == 'Pathways':
                    a = item['Gene-commentary_comment']
                    for b in a:
                        c = b['Gene-commentary_text']
                        pathwaylist.append(c)
        return pathwaylist
    except:
        return None

def PMID(record):
    pmidlist = []
    try:
        pmid = record['Entrezgene_comments']
        for item in pmid:
            for key,value in item.iteritems():
                if key == 'Gene-commentary_refs':
                    for item in value:
                        for a,b in item.iteritems():
                            if a == 'Pub_pmid':
                                c = b['PubMedId']
                                pmidlist.append(c)
        return pmidlist
    except:
        return None



def Interaction(record):
    interactionlist = []
    try:
        interaction = record['Entrezgene_comments']
        for item in interaction:
            for key,value in item.iteritems():
                if key == 'Gene-commentary_comment':
                    if item['Gene-commentary_heading'] == 'Interactions':
                        a = item['Gene-commentary_comment']
                        for b in a:
                            c = b['Gene-commentary_comment']
                            for d in c:
                                e = d['Gene-commentary_source']
                                for f in e:
                                    for k,v in f.iteritems():
                                        if k == 'Other-source_anchor':
                                            interactionlist.append(v)
        return interactionlist                                                        
    except:
        return None 

def Conserved(record):
    conservedlist = []
    try:
        conserved = record['Entrezgene_comments']
        for item in conserved:
            for key,value in item.iteritems():
                if key == 'Gene-commentary_comment':
                    for item in value:
                        for a,b in item.iteritems():
                            if a == 'Gene-commentary_products':
                                for c in b:
                                    for aa,bb in c.iteritems():
                                        if aa == 'Gene-commentary_products':                                 
                                            d = c['Gene-commentary_products']
                                            for e in d:
                                                f = e['Gene-commentary_comment']
                                                for g in f:
                                                    if g['Gene-commentary_heading'] == 'Conserved Domains':
                                                        try:
                                                            h = g['Gene-commentary_comment'],'\n'
                                                            for i in h:
                                                                for j in i:
                                                                    try:
                                                                        k = j['Gene-commentary_source']
                                                                        for l in k:
                                                                            try:
                                                                                j = l['Other-source_anchor']
                                                                                conservedlist.append(j)
                                                                            except:
                                                                                continue
                                                                    except:
                                                                        continue
                                                        except:
                                                            continue
        return conservedlist
    except:
        print 'Whoops!','\n'
        pass


def Function(record):
    functionlist = []
    try:
        function = record['Entrezgene_properties']
        for item in function:
            for key,value in item.iteritems():
                if key == 'Gene-commentary_comment':
                    for item in value:
                        for key,value in item.iteritems():
                            if value == 'Function':
                                a = item['Gene-commentary_comment']
                                for b in a:
                                    c = b['Gene-commentary_source']
                                    for d in c:
                                        e = d['Other-source_anchor']
                                        functionlist.append(e)
        return functionlist
    except:
        return None


def Component(record):
    componentlist = []
    try:
        component = record['Entrezgene_properties']
        for item in component:
            for key,value in item.iteritems():
                if key == 'Gene-commentary_comment':
                    for item in value:
                        for key,value in item.iteritems():
                            if value == 'Component':
                                a = item['Gene-commentary_comment']
                                for b in a:
                                    c = b['Gene-commentary_source']
                                    for d in c:
                                        e = d['Other-source_anchor']
                                        componentlist.append(e)
        return componentlist
    except:
        return None

def Process(record):
    processlist = []
    try:
        process = record['Entrezgene_properties']
        for item in process:
            for key,value in item.iteritems():
                if key == 'Gene-commentary_comment':
                    for item in value:
                        for key,value in item.iteritems():
                            if value == 'Process':
                                a = item['Gene-commentary_comment']
                                for b in a:
                                    c = b['Gene-commentary_source']
                                    for d in c:
                                        e = d['Other-source_anchor']
                                        processlist.append(e)
        return processlist
    except:
        return None
    



if __name__ == '__main__':
    main()



