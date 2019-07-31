import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import dash_table
import pandas as pd
import os
import ursgal

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

protein_dict = {}
db_fasta = os.path.join('databases', 'Haloferax_volcanii_ArcPP_20190606_uniprot.fasta')
with open(db_fasta, 'r') as db_input:
    for fasta_id, db_sequence in ursgal.ucore.parse_fasta(db_input):
        protein_dict[fasta_id] = db_sequence

coverage_dict = {}
df = pd.read_csv(os.path.join('results', 'ArcPP_results.csv'))
for index, row in df.iterrows():
    prot = row['Protein ID']
    try:
        start = int(row['Sequence Start'])
        stop = int(row['Sequence Stop'])
    except:
        continue
    if prot not in coverage_dict.keys():
        coverage_dict[prot] = set()
    for aa in range(start, stop+1):
        coverage_dict[prot].add(aa)

app.layout = html.Div([
    dash_table.DataTable(
        id='datatable-interactivity',
        columns=[
            {"name": i, "id": i, "deletable": True} for i in df.columns
        ],
        data=df.to_dict('records'),
        editable=True,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_selectable="single",
        row_deletable=False,
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 20,
    ),
    # html.Div(id='datatable-interactivity-container')
])

# app.layout = html.Div([
#     dcc.Input(id='my-id', value='initial value', type='text'),
#     html.Div(id='my-div')
# ])

# @app.callback(
#     Output('datatable-interactivity-container', "children"),
#     [Input('datatable-interactivity', "derived_virtual_selected_rows")]
# )

# # @app.callback(
# #     Output(component_id='my-div', component_property='children'),
# #     [Input(component_id='my-id', component_property='value')]
# # )
# def update_output_div(derived_virtual_selected_rows):
#     if derived_virtual_selected_rows is None:
#         derived_virtual_selected_rows = []
        
#     if derived_virtual_selected_rows == []:
#         protein = 'No protein selected'
#         sequence_list = ['']

#     else:
#         protein = df.iloc[derived_virtual_selected_rows]['Protein ID'].values[0]
#         sequence_list = []

#         s = ''
#         for n, aa in enumerate(protein_dict[protein]):
#             if n in coverage_dict[protein]:
#                 s += '*{0}*'.format(aa)
#             else:
#                 s += aa
#             if (n + 1) % 80 == 0:
#                 sequence_list.append(s)
#                 s = ''

#     return '''

#         >{0}
#         {1}
#         '''.format(protein, '\n'.join(sequence_list))




if __name__ == '__main__':
    app.run_server(debug=True)
