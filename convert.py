import pdftables_api

c = pdftables_api.Client('26xisrk36qt7')
c.xlsx('ALL.pdf', 'ALL.xlsx')
c.xlsx('AML.pdf', 'AML.xlsx')
