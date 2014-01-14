import html
import datetime

DEFAULT_REPORT_FILE = 'annokey_report.html'
NCBI_GENE_ENTRY_URL = "http://www.ncbi.nlm.nih.gov/gene/"
GENE_RIF_URL = "http://www.ncbi.nlm.nih.gov/gene?db=gene&report=generif&term="
PUBMED_URL = "http://www.ncbi.nlm.nih.gov/gene/pubmed?LinkName=gene_pubmed&form_uid="

style_css = """
body
{
font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
width:800px;
margin:0px auto 0px auto;
text-align:justify
}
#hitstable
{
font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
width:100%;
border-collapse:collapse;
}
#hitstable td, #hitstable th 
{
font-size:1em;
border:1px solid #98bf21;
padding:3px 7px 2px 7px;
text-align:left;
}
#hitstable th 
{
font-size:1.1em;
text-align:left;
padding-top:5px;
padding-bottom:4px;
background-color:#A7C942;
color:#ffffff;
}
#hitstable tr.alt td 
{
color:#000000;
background-color:#EAF2D3;
}
"""

def init_report_page(command_line_text):
    page = html.HTML('html')
    title = 'Annokey Search Report'
    page.title(title)
    page.h1(title)
    page.style(style_css)
    report_command_date(page, command_line_text)
    return page

def report_command_date(page, command_line_text):
    list = page.ul
    command_item = list.li
    command_item.b("Command: ")
    command_item.text(command_line_text)
    date_item = list.li
    date_item.b("Date: ")
    today = datetime.date.today()
    date_item.text(today.strftime('%d %b %Y'))

def report_hits(gene_name, hits, report_page):
    # don't report a gene if it has no hits
    if len(hits) > 0:
        report_page.h2(gene_name)
        ncbi_record_id = str(hits[0].database_record_id)
        para = report_page.p
        ncbi_gene_url = NCBI_GENE_ENTRY_URL + ncbi_record_id
        gene_rif_url = GENE_RIF_URL + ncbi_record_id
        pubmed_url = PUBMED_URL + ncbi_record_id
        list = para.ul
        item = list.li
        item.a("NCBI Gene", href=ncbi_gene_url)
        item = list.li
        item.a("GeneRIF", href=gene_rif_url)
        item = list.li
        item.a("Pubmed", href=pubmed_url)
        table = report_page.table(id="hitstable")
        row = table.tr
        row.th("Search Term")
        row.th("Rank")
        row.th("Hits")
        alternate_row = False
        for h in hits:
            if alternate_row:
                row = table.tr(klass = "alt")
                alternate_row = False 
            else:
                row = table.tr
                alternate_row = True
            row.td(h.search_term)
            row.td(str(h.rank))
            td = row.td('; '.join(h.fields))

def write_report(report_filename, report_page):
    with open(report_filename, "w") as report_file:
        report_file.write(str(report_page))
