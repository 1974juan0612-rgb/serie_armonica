import os
from fpdf import FPDF

class PDF(FPDF):
    def header(self):
        pass
    def footer(self):
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)
        self.cell(0, 10, f'{self.page_no()}/{{nb}}', 0, 0, 'C')
    def chapter_title(self, title):
        self.set_font('Helvetica', 'B', 16)
        self.multi_cell(0, 8, title)
        self.ln(4)
    def section_title(self, title):
        self.set_font('Helvetica', 'B', 12)
        self.multi_cell(0, 7, title)
        self.ln(2)
    def body_text(self, text):
        self.set_font('Helvetica', '', 10)
        self.multi_cell(0, 5.5, text)
        self.ln(2)
    def math_block(self, text):
        self.set_font('Courier', '', 9)
        self.multi_cell(0, 5, text)
        self.ln(2)

OUT = r'C:\Users\famil\Desktop\serie_armonica\documentos'
SRC = OUT

articles = [
    ('articulo_escala_geometrica.md', 'escala_geometrica.pdf', 'La escala geometrica del colapso en redes FCC con DNLS'),
    ('articulo_modelo_unificado.md', 'modelo_unificado.pdf', 'Modelo unificado de red FCC con cero parametros libres'),
    ('articulo_barrera_cuantico_clasico.md', 'barrera_cuantico_clasico.pdf', 'Emergencia de la barrera cuantico-clasico'),
]

for src_name, pdf_name, short_title in articles:
    path = os.path.join(SRC, src_name)
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    pdf = PDF()
    pdf.alias_nb_pages()
    pdf.add_page()

    # Title
    pdf.set_font('Helvetica', 'B', 18)
    title = lines[0].replace('#', '').strip() if lines[0].startswith('#') else short_title
    pdf.multi_cell(0, 10, title)
    pdf.ln(3)

    # Author/date
    for line in lines[1:5]:
        if line.startswith('**') or line.startswith('*'):
            pdf.set_font('Helvetica', 'I', 10)
            pdf.multi_cell(0, 6, line.strip('* \n\r'))
            pdf.ln(1)

    pdf.ln(5)

    # Body
    i = 5
    while i < len(lines):
        line = lines[i].rstrip()
        i += 1

        if not line:
            continue

        # Section headers
        if line.startswith('## '):
            pdf.section_title(line.replace('## ', '').strip())
            continue
        if line.startswith('### '):
            pdf.set_font('Helvetica', 'B', 11)
            pdf.multi_cell(0, 6, line.replace('### ', '').strip())
            pdf.ln(1)
            continue

        # Skip standalone # title line (already handled)
        if line.startswith('# ') and line.strip() == title:
            continue

        # Render as body
        pdf.body_text(line)

    pdf.output(os.path.join(OUT, pdf_name))
    print(f'Creado: {OUT}\\{pdf_name}')

print('Listo. 3 PDFs generados.')
