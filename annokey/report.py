'''
Generate a HTML report of the results of an
Annokey search.

The report file contains a much more detailed
presentation of the search results than what
is shown in the annotation on the CSV file.
'''

import html
import datetime

DEFAULT_REPORT_FILE = 'annokey_report.html'
NCBI_GENE_ENTRY_URL = "http://www.ncbi.nlm.nih.gov/gene/"
GENE_RIF_URL = "http://www.ncbi.nlm.nih.gov/gene?db=gene&report=generif&term="
PUBMED_URL = \
    "http://www.ncbi.nlm.nih.gov/gene/pubmed?LinkName=gene_pubmed&from_uid="

STYLE_CSS = """
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
.detail{
   display:none;
}
.table_col_term{
   width: 45%;
}
.table_col_rank{
   width: 10%;
}
.table_col_fields{
   width: 45%;
}
"""

JAVASCRIPT = """
    function toggle_visibility(id) {
       var e = document.getElementById(id);
       if(e.style.display == 'block')
          e.style.display = 'none';
       else
          e.style.display = 'block';
    }
"""

ANNOKEY_URL = 'http://bjpop.github.io/annokey/'


def init_report_page(hostname, working_directory, command_line_text):
    '''Generare the start of the output page, including heading
    and initial part of the body.
    '''
    page = html.HTML('html')
    head = page.head
    title = 'Annokey Search Report'
    head.title(title)
    body = page.body
    body.h1(title)
    body.style(STYLE_CSS)
    body.script(JAVASCRIPT)
    report_meta_data(body, hostname, working_directory, command_line_text)
    return page

def report_meta_data(page, hostname, working_directory, command_line_text):
    '''Generate the part of the page containing meta information about
    the annokey search, such as the command that was run, the
    hostname of the machine on which it was run, the date, a link to
    the annokey homepage.
    '''
    list_element = page.ul
    annokey_item = list_element.li
    annokey_item.b("annokey homepage: ")
    annokey_item.a(ANNOKEY_URL, href=ANNOKEY_URL)
    hostname_item = list_element.li
    hostname_item.b("hostname: ")
    hostname_item.text(hostname)
    directory_item = list_element.li
    directory_item.b("directory: ")
    directory_item.code(working_directory)
    command_item = list_element.li
    command_item.b("command: ")
    command_item.code(command_line_text)
    date_item = list_element.li
    date_item.b("date: ")
    today = datetime.date.today()
    date_item.text(today.strftime('%d %b %Y'))

def report_hits(gene_count, gene_name, gene_db_id, hits, report_page):
    '''Generate the HTML for all the search hits for all
    the genes with at least one hit.
    '''
    # don't report a gene if it has no hits
    if hits:
        report_page.h2(gene_name)
        ncbi_gene_url = NCBI_GENE_ENTRY_URL + gene_db_id
        gene_rif_url = GENE_RIF_URL + gene_db_id
        pubmed_url = PUBMED_URL + gene_db_id
        list_element = report_page.ul
        item = list_element.li
        item.a("NCBI Gene", href=ncbi_gene_url)
        item = list_element.li
        item.a("GeneRIF", href=gene_rif_url)
        item = list_element.li
        item.a("Pubmed", href=pubmed_url)
        sorted_hits = sorted(hits)
        make_hit_table(report_page, hits, sorted_hits)
        make_detailed_match_lists(gene_count, report_page, hits,
            sorted_hits, gene_name)


def make_hit_table(report_page, hits, sorted_hits):
    '''Generate a table which summarises the key term hits
    for a given gene.
    '''
    table = report_page.table(id="hitstable")
    row = table.tr
    row.th("Search Term", klass="table_col_term")
    row.th("Rank", klass="table_col_rank")
    row.th("Fields", klass="table_col_fields")
    alternate_row = False
    for (rank, term) in sorted_hits:
        fields = hits[(rank, term)]
        if alternate_row:
            row = table.tr(klass="alt")
            alternate_row = False
        else:
            row = table.tr
            alternate_row = True
        row.td(term)
        row.td(str(rank))
        field_info = []
        for field in fields:
            field_count = len(fields[field])
            field_info.append("{}({})".format(field, field_count))
        row.td('; '.join(field_info))


def make_detailed_match_lists(gene_count, report_page, hits, sorted_hits,
    gene_name):
    '''For a given gene, for each key term, generate a detailed
    report of each occurence of the term in the corresponding
    matched field.
    '''
    _ = report_page.br
    gene_div_id = "gene" + str(gene_count)
    report_page.button("hide/show details for {}".format(gene_name),
        onclick="toggle_visibility('{}');".format(gene_div_id))
    div = report_page.div(id=gene_div_id, klass='detail')
    for (rank, term) in sorted_hits:
        fields = hits[(rank, term)]
        div.h3(term)
        for field, matches in fields.items():
            make_match_list(div, field, matches)
    _ = div.hr # underscore makes pylint happy

def make_match_list(div, field, matches):
    '''Generate a list of matches and highlight
    the matched key term in the output.
    '''
    div.h4(field)
    with div.ol as list_element:
        for match in matches:
            with list_element.li as list_element_item:
                start = 0
                context = to_string(match.context)
                for (span_low, span_high) in match.spans:
                    before_span = context[start:span_low]
                    list_element_item.text(before_span)
                    list_element_item.b.u(
                        context[span_low:span_high])
                    start = span_high
                list_element_item.text(context[start:])


def to_string(text):
    '''If a piece of text is unicode, turn it into an ASCII string.'''
    if isinstance(text, unicode):
        return text.encode('utf8')
    else:
        return text

def write_report(report_filename, report_page):
    '''Save the HTML report file.'''
    with open(report_filename, "w") as report_file:
        report_file.write(str(report_page))
