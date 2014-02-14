import html
import datetime

DEFAULT_REPORT_FILE = 'annokey_report.html'
NCBI_GENE_ENTRY_URL = "http://www.ncbi.nlm.nih.gov/gene/"
GENE_RIF_URL = "http://www.ncbi.nlm.nih.gov/gene?db=gene&report=generif&term="
PUBMED_URL = "http://www.ncbi.nlm.nih.gov/gene/pubmed?LinkName=gene_pubmed&from_uid="

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

javascript = """
    function toggle_visibility(id) {
       var e = document.getElementById(id);
       if(e.style.display == 'block')
          e.style.display = 'none';
       else
          e.style.display = 'block';
    }
"""

# we need a unique ID for each gene in the entire output file
# we use it to generate a unique div id, which enables us
# to hide and show parts of the page
gene_count = 0

def init_report_page(hostname, working_directory, command_line_text):
    page = html.HTML('html')
    head = page.head
    title = 'Annokey Search Report'
    head.title(title)
    body = page.body
    body.h1(title)
    body.style(style_css)
    body.script(javascript)
    report_meta_data(body, hostname, working_directory, command_line_text)
    return page

def report_meta_data(page, hostname, working_directory, command_line_text):
    list = page.ul
    annokey_item = list.li
    annokey_item.b("annokey homepage: ")
    annokey_item.a("http://bjpop.github.io/annokey/", href="http://bjpop.github.io/annokey/")
    hostname_item = list.li
    hostname_item.b("hostname: ")
    hostname_item.text(hostname)
    directory_item = list.li
    directory_item.b("directory: ")
    directory_item.code(working_directory)
    command_item = list.li
    command_item.b("command: ")
    command_item.code(command_line_text)
    date_item = list.li
    date_item.b("date: ")
    today = datetime.date.today()
    date_item.text(today.strftime('%d %b %Y'))

def report_hits(gene_name, gene_db_id, hits, report_page):
    # don't report a gene if it has no hits
    global gene_count
    if hits:
        report_page.h2(gene_name)
        ncbi_gene_url = NCBI_GENE_ENTRY_URL + gene_db_id
        gene_rif_url = GENE_RIF_URL + gene_db_id
        pubmed_url = PUBMED_URL + gene_db_id
        list = report_page.ul
        item = list.li
        item.a("NCBI Gene", href=ncbi_gene_url)
        item = list.li
        item.a("GeneRIF", href=gene_rif_url)
        item = list.li
        item.a("Pubmed", href=pubmed_url)
        table = report_page.table(id="hitstable")
        row = table.tr
        row.th("Search Term", klass="table_col_term")
        row.th("Rank", klass="table_col_rank")
        row.th("Fields", klass="table_col_fields")
        alternate_row = False
        sorted_hits = sorted(hits)
        for (rank, term) in sorted_hits:
            fields = hits[(rank, term)]
            if alternate_row:
                row = table.tr(klass = "alt")
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
        report_page.br
        gene_div_id = "gene" + str(gene_count)
        gene_count += 1
        report_page.button("hide/show details for {}".format(gene_name), onclick="toggle_visibility('{}');".format(gene_div_id))
        div = report_page.div(id=gene_div_id, klass='detail')
        for (rank, term) in sorted_hits:
            fields = hits[(rank, term)]
            div.h3(term)
            for field in fields:
                div.h4(field)
                with div.ol as list:
                    for (before, match, after) in fields[field]:
                        with list.li as list_item:
                            list_item.text(to_string(before))
                            list_item.b(to_string(match))
                            list_item.text(to_string(after))
        div.hr

def to_string(text):
    if isinstance(text, unicode):
        return text.encode('utf8')
    else:
        return text
                    

def write_report(report_filename, report_page):
    with open(report_filename, "w") as report_file:
        report_file.write(str(report_page))
